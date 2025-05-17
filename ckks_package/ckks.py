from sage.all import (
    ceil,
    divisors,
    factorial,
    floor,
    imag,
    Integer,
    log,
    real,
    round,
)
import numpy as np
from .bit_rev import bit_rev_vector
from .fast_dft import get_grouped_E, get_grouped_F, Multidiags
from .lattice_estimator.estimator import LWE, ND
from .poly import Poly


class CKKS:
    """
    Class for CKKS scheme.
    """

    @classmethod
    def config(cls, N, n, L_boot, q0, p, delta, print_messages=False):
        """
        Configure basic CKKS parameters.

        Args:
            N (int):
                Ring degree (must be a power of two).
            n (int):
                Number of slots (must be a strict divisor of N).
            L_boot (int):
                Maximal level during bootstrapping.
            q0 (int):
                Base modulus.
            p (int):
                Scaling factor outside of bootstrapping.
            delta (int):
                Scaling factor during bootstrapping.
            print_messages (bool, optional):
                Whether to print information messages. Defaults to False.

        Raises:
            ValueError:
                If N is not a power of two.
            ValueError:
                If n is not a strict divisor of N.
            ValueError:
                If p and delta are not powers of the same base.
        """
        # Generating matrices required for encoding and decoding
        get_grouped_E(n, 1, inverse=True)
        get_grouped_E(n, 1, inverse=False)

        if N & N - 1 != 0:
            raise ValueError("N must be a power of two.")
        if N % n != 0 or n >= N:
            raise ValueError("n must be a strict divisor of N.")

        are_same_powers = False
        for d in divisors(p):
            if isinstance(log(p, d), (int, Integer)) and isinstance(
                log(delta, d), (int, Integer)
            ):
                are_same_powers = True
                break
        if not are_same_powers:
            raise ValueError("p and delta must be powers of the same base.")

        cls.N = N  # Ring degree
        cls.n = n  # Number of slots
        cls.q0 = q0  # Base modulus
        cls.p = p  # Scaling factor outside of bootstrapping
        cls.delta = delta  # Scaling factor during bootstrapping

        cls.L = floor(
            L_boot * log(delta, p) + 1
        )  # Maximal level outside of bootstrapping
        cls.L_boot = L_boot  # Maximal level for bootstrapping

        cls.moduli = [
            q0 * p**i for i in range(cls.L + 1)
        ]  # The ct moduli used outside of bootstrapping
        cls.moduli_boot = [
            p * q0 * delta**i for i in range(L_boot + 1)
        ]  # The ct moduli used during bootstrapping

        if print_messages:
            print("The CKKS configuration is done!")

    # Key generation

    @classmethod
    def key_gen(
        cls, h=None, sk=None, q=None, P=None, sigma=None, print_messages=False
    ):
        """
        Generate the secret key, public key and evaluation key for the scheme.

        Args:
            h (int, optional):
                Hamming weight of the secret key. Defaults to
                2 ** min(6, (log(N, 2) // 2)).
            sk (Poly, optional):
                Secret key for scheme. Defaults to a ternary polynomial of
                Hamming weight equal to h.
            q (int, optional):
                Largest modulus for the evaluation key. Defaults to the largest
                modulus supported by the scheme.
            P (int, optional):
                Extra factor added to the evaluation key modulus. Defaults to
                the largest modulus accepted by the scheme.
            sigma (float, optional):
                Standard deviation for error polynomials. Defaults to 3.2.
            print_messages (bool, optional):
                Whether to print information messages. Defaults to False.

        Raises:
            RuntimeError:
                If CKKS is not configured before key generation.
        """
        try:
            cls.N
        except AttributeError:
            raise RuntimeError("Ask Alice to configure CKKS first.")

        if sk is not None:
            cls.h = sum(sk.coeffs[i] != 0 for i in range(cls.N))
        elif h is not None:
            cls.h = h
        else:
            cls.h = 2 ** min(6, 2 ** (log(cls.N, 2) // 2))

        cls.P = cls.moduli_boot[-1] if P is None else P
        cls.sigma = 3.2 if sigma is None else sigma

        q_sk = cls.P * cls.moduli_boot[-1]  # Modulus for secret key
        if sk is not None:
            cls.sk = sk % q_sk
        else:
            cls.sk = Poly.get_random_ternary1(cls.N, q_sk, cls.h)

        e = Poly.get_random_normal(cls.N, cls.moduli_boot[-1], cls.sigma)
        pk1 = Poly.get_random_uniform(cls.N, cls.moduli_boot[-1])
        pk0 = e - pk1 * cls.sk
        cls.pk = (pk0, pk1)

        q_evk = cls.P * (
            cls.moduli_boot[-1] if q is None else q
        )  # Modulus for evk
        e = Poly.get_random_normal(cls.N, q_evk, cls.sigma)
        evk1 = Poly.get_random_uniform(cls.N, q_evk)
        evk0 = e + cls.P * cls.sk**2 - evk1 * cls.sk
        cls.evk = (evk0, evk1)

        cls.galois_swk_dict = {}  # Dictionary containing the swk for Galois

        if print_messages:
            print("The key generation is done!")

    @classmethod
    def _check_key_gen(cls):
        """
        Check if CKKS is configured and keys have been generated.

        Raises:
            RuntimeError:
                If CKKS is not yet configured.
            RuntimeError:
                If keys have not yet been generated.
        """
        try:
            cls.N
        except:
            raise RuntimeError(
                "Ask Alice to configure CKKS and to generate the keys first."
            )
        try:
            cls.evk
        except:
            raise RuntimeError("Ask Alice to generate the keys first.")

    @classmethod
    def get_swk(cls, sk_x, sk, q=None, P=None):
        """
        Generate a switching key to switch from secret key sk_x to sk.

        Args:
            sk_x (Poly):
                Secret key to switch from.
            sk (Poly):
                Secret key to switch to.
            q (int, optional):
                Largest modulus for the switching key. Defaults to the largest
                modulus supported by the scheme.
            P (int, optional):
                Extra factor added to the modulus. Defaults to cls.P.

        Returns:
            tuple:
                Switching key (swk0, swk1).
        """
        cls._check_key_gen()

        if q is None:
            q = cls.moduli_boot[-1]
        if P is None:
            P = cls.P

        q_evk = q * P
        e = Poly.get_random_normal(cls.N, q_evk, cls.sigma)
        swk1 = Poly.get_random_uniform(cls.N, q_evk)  # Slow
        swk0 = e + P * sk_x - swk1 * sk
        return (swk0, swk1)

    @classmethod
    def get_galois_swk(cls, k, sk=None, q=None, P=None):
        """
        Get switching key for applying Galois automorphism.

        Args:
            k (int):
                Exponent of the Galois automorphism.
            sk (Poly, optional):
                Secret key. Required if the switching key is not already
                generated.
            q (int, optional):
                Largest modulus for the switching key.
            P (int, optional):
                Extra factor added to the modulus.

        Returns:
            tuple:
                Switching key for Galois automorphism.

        Raises:
            RuntimeError:
                If sk is None and the switching key is not already generated.
        """
        k = k % (2 * cls.N)
        cls._check_key_gen()
        if k in cls.galois_swk_dict:
            return cls.galois_swk_dict[k]
        if sk is None:
            raise RuntimeError("Secret key sk required.")
        sk_x = sk.galois(k)
        swk = cls.get_swk(sk_x, sk, q, P)
        cls.galois_swk_dict[k] = swk
        return swk

    # Security estimation

    @classmethod
    def get_security(
        cls,
        N=None,
        q=None,
        sigma=None,
        sk_plus=None,
        sk_minus=None,
        check_primal_hybrid=False,
    ):
        """
        Use the LWE estimator to compute the security level.

        Args:
            N (int, optional):
                Ring degree. Defaults to cls.N.
            q (int, optional):
                Largest modulus used for switching keys. Defaults to the
                largest modulus supported by the scheme, multiplied by P.
            sigma (float, optional):
                Standard deviation for error polynomials. Defaults to
                cls.sigma.
            sk_plus (int, optional):
                Number of coefficients equal to 1 in the secret key. Defaults
                to the number of coefficients equal to 1 in cls.sk.
            sk_minus (int, optional):
                Number of coefficients equal to -1 in the secret key. Defaults
                to the number of coefficients equal to -1 in cls.sk.
            check_primal_hybrid (bool, optional):
                Whether to check the primal hybrid security level. Defaults
                to False.

        Returns:
            float:
                Logarithm (base 2) of estimated number of operations required
                to break the scheme.
        """
        if None in (N, q, sigma, sk_plus, sk_minus):
            cls._check_key_gen()

        if N is None:
            N = cls.N
        if q is None:
            q = cls.moduli_boot[-1] * cls.P
        if sigma is None:
            sigma = cls.sigma
        if sk_plus is None or sk_minus is None:
            sk_plus = 0
            sk_minus = 0
            for i in range(cls.N):
                sk_plus += cls.sk.coeffs[i] == 1
                sk_minus += cls.sk.coeffs[i] == -1

        params = LWE.Parameters(
            N,
            q,
            Xs=ND.SparseTernary(n=N, p=sk_plus, m=sk_minus),
            Xe=ND.DiscreteGaussian(n=N, stddev=sigma),
        )

        s0 = LWE.dual(params)["rop"]
        s1 = LWE.dual_hybrid(params)["rop"]
        s2 = LWE.primal_usvp(params)["rop"]
        if check_primal_hybrid:
            s3 = LWE.primal_hybrid(params)["rop"]
        else:
            s3 = s0

        return round(log(min(s0, s1, s2, s3), 2), 2)

    # Initialization of ciphertexts

    @classmethod
    def _check_modulus(cls, q, is_boot=False):
        """
        Check if the current parameters support the given modulus q.

        Args:
            q (int):
                Modulus to check.
            is_boot (bool, optional):
                Whether the given modulus q is used during bootstrapping.
                Defaults to False.

        Raises:
            ValueError:
                If modulus q is not supported by the current CKKS parameters.
        """
        if q == 0:
            return
        if (is_boot and q not in cls.moduli_boot) or (
            not is_boot and q not in cls.moduli
        ):
            raise ValueError(
                "Modulus q not supported by current CKKS parameters."
            )

    def __init__(self, b, a, is_boot=False):
        """
        Initialize a CKKS ciphertext with tuple (b, a).

        Args:
            b (Poly):
                Polynomial b.
            a (Poly):
                Polynomial a.
            is_boot (bool, optional):
                Whether the ciphertext is used during bootstrapping.
                Defaults to False.

        Raises:
            ValueError:
                If b or a have incorrect degree N.
            ValueError:
                If b and a have different moduli q.
        """
        if b.N != self.N or a.N != self.N:
            raise ValueError("b or a has wrong degree N.")
        if b.q != a.q:
            raise ValueError("b and a have different moduli q.")

        self.is_boot = is_boot

        if a.q == 0:
            q = self.moduli_boot[-1] if is_boot else self.moduli[-1]
            self.b, self.a = b % q, a % q
            self.q = q
            self.l = self.L_boot if is_boot else self.L
        else:
            self._check_modulus(a.q, is_boot)
            self.b, self.a = b, a
            self.q = a.q
            self.l = (
                self.moduli_boot.index(a.q)
                if is_boot
                else self.moduli.index(a.q)
            )

    # Encoding and decoding

    @classmethod
    def encode(cls, z, is_boot=False):
        """
        Encode a complex vector z into an integer polynomial. This process
        requires applying the inverse of a variant of the DFT matrix, and
        scaling the result by the factor p or delta to achieve the necessary
        precision.

        Args:
            z (np.ndarray or np.complex128):
                A complex vector (of length 1 if z is an np.complex128) to
                encode.
            is_boot (bool, optional):
                Whether the polynomial is used during bootstrapping. Defaults
                to False.

        Returns:
            Poly:
                Integer polynomial encoding the vector z.

        Raises:
            ValueError:
                If the length of z is not a strict divisor of N.
        """
        try:
            n = len(z)
            if n == 1:
                z = z[0]
        except TypeError:
            n = 1

        if cls.N % n != 0 or n >= cls.N:
            raise ValueError("n must be a strict divisor of N.")

        scaling_factor = cls.delta if is_boot else cls.p
        q = cls.moduli_boot[-1] if is_boot else cls.moduli[-1]

        if n == 1:
            # Single value, representing a vector whose entries are all equal
            # to this value.
            a = int(round(scaling_factor * real(z)))
            b = int(round(scaling_factor * imag(z)))
            return Poly.get_constant(a, cls.N, q) + b * Poly.get_monomial(
                cls.N // 2, cls.N, q
            )

        grouped_iE = get_grouped_E(n, 1, inverse=True)

        w = bit_rev_vector(z, num_bits=log(n, 2))
        for A in grouped_iE:
            w = A.BSGS_mult(w)
        w *= scaling_factor

        pt0_temp = np.round(real(w)).astype(int)
        pt1_temp = np.round(imag(w)).astype(int)
        pt0, pt1 = [np.zeros(cls.N // 2).astype(int) for _ in range(2)]

        l = cls.N // (2 * n)
        for i in range(n):
            pt0[l * i] = pt0_temp[i]
            pt1[l * i] = pt1_temp[i]

        return Poly(np.concatenate((pt0, pt1)), cls.N, q)

    @classmethod
    def decode(cls, pt, n=None, is_boot=False):
        """
        Decode an integer polynomial to a complex vector of length n. This
        process is effectively the inverse of the previously described encode
        method.

        Args:
            pt (Poly):
                Integer polynomial to decode.
            n (int, optional):
                Length of the complex vector. Defaults to cls.n.
            is_boot (bool, optional):
                Whether the input polynomial was used during bootstrapping.
                Defaults to False.

        Returns:
            np.ndarray:
                Complex vector z of length n.

        Raises:
            ValueError:
                If n is not a strict divisor of N.
        """
        if n is None:
            n = cls.n
        elif cls.N % n != 0 or n >= cls.N:
            raise ValueError("n must be a strict divisor of N.")

        grouped_E = get_grouped_E(n, 1, inverse=False)

        l = cls.N // (2 * n)
        coeffs = pt.get_symmetric_coeffs()
        pt0 = np.array([int(coeffs[l * i]) for i in range(n)])
        pt1 = np.array([int(coeffs[l * i]) for i in range(n, 2 * n)])

        w = pt0.astype(np.complex128) + 1j * pt1.astype(np.complex128)
        for A in grouped_E:
            w = A.BSGS_mult(w)
        scaling_factor = cls.delta if is_boot else cls.p
        w /= scaling_factor

        return bit_rev_vector(w, log(n, 2))

    # Encryption and decryption

    @classmethod
    def enc_poly_with_sk(cls, pt, sk, is_boot=False):
        """
        Encrypt a plaintext polynomial with the secret key.

        Args:
            pt (Poly):
                Plaintext polynomial to encrypt.
            sk (Poly):
                Secret key.
            is_boot (bool, optional):
                Whether the output ciphertext should be used during
                bootstrapping. Defaults to False.

        Returns:
            CKKS:
                Encrypted ciphertext.
        """
        cls._check_modulus(pt.q, is_boot)
        if pt.q != 0:
            q = pt.q
        elif is_boot == False:
            q = cls.moduli[-1]
        else:
            q = cls.moduli_boot[-1]
        e = Poly.get_random_normal(cls.N, q, cls.sigma)
        a = Poly.get_random_uniform(cls.N, q)
        b = pt + e - a * sk
        return cls(b, a, is_boot)

    @classmethod
    def enc_poly_with_pk(cls, pt, pk):
        """
        Encrypt a plaintext polynomial with the public key.

        Args:
            pt (Poly):
                Plaintext polynomial to encrypt.
            pk (tuple):
                Public key (pk0, pk1).

        Returns:
            CKKS:
                Encrypted ciphertext.
        """
        cls._check_modulus(pt.q)
        q = cls.moduli[-1] if pt.q == 0 else pt.q
        v = Poly.get_random_ternary2(cls.N, q)
        e0, e1 = [
            Poly.get_random_normal(cls.N, q, cls.sigma) for _ in range(2)
        ]
        return cls(pt + v * pk[0] + e0, v * pk[1] + e1, is_boot=False)

    @classmethod
    def enc_poly_without_error(cls, pt, is_boot=False):
        """
        Encrypt a plaintext polynomial without error (insecure).

        Args:
            pt (Poly or int):
                Plaintext polynomial or integer to encrypt.
            is_boot (bool, optional):
                Whether the ciphertext is used during bootstrapping.
                Defaults to False.

        Returns:
            CKKS:
                Encrypted ciphertext.
        """
        q = cls.moduli_boot[-1] if is_boot else cls.moduli[-1]
        if isinstance(pt, (int, Integer)):
            a = Poly.get_constant(0, cls.N, q)
            b = Poly.get_constant(pt, cls.N, q)
            return cls(b, a, is_boot)
        cls._check_modulus(pt.q, is_boot)
        q = q if pt.q == 0 else pt.q
        a = Poly.get_constant(0, cls.N, q)
        b = pt % q
        return cls(b, a, is_boot)

    def dec_to_poly(self, sk):
        """
        Decrypt the ciphertext to a plaintext polynomial.

        Args:
            sk (Poly):
                Secret key.

        Returns:
            Poly:
                Decrypted plaintext polynomial.
        """
        return self.b + self.a * sk

    # Modular reduction, lifting and rescaling

    def __mod__(self, q):
        """
        Compute self % q.

        Args:
            q (int):
                Modulus.

        Returns:
            CKKS:
                Resulting ciphertext.
        """
        return self.__class__(self.b % q, self.a % q, self.is_boot)

    def __imod__(self, q):
        """
        In place modulo operation with modulus q.

        Args:
            q (int):
                Modulus.

        Returns:
            CKKS:
                Resulting ciphertext.
        """
        return self % q

    def lift(self, q=None):
        """
        Lift the ciphertext to a higher modulus.

        Args:
            q (int, optional):
                The new modulus to lift the ciphertext to. Defaults to the
                largest modulus supported by the CKKS scheme.

        Returns:
            CKKS:
                Resulting ciphertext.
        """
        if q is None and self.is_boot:
            q = self.moduli_boot[-1]
        elif q is None:
            q = self.moduli[-1]
        return __class__(self.b.lift(q), self.a.lift(q), self.is_boot)

    def rescale(self, l=1):
        """
        Rescale by dividing by p**l or delta**l.

        Args:
            l (int, optional):
                Number of levels to rescale. Defaults to 1.

        Returns:
            CKKS:
                Rescaled ciphertext.

        Raises:
            ValueError:
                If attempting to rescale to a negative level.
        """
        if self.l < l:
            raise ValueError("Cannot rescale to a negative level.")
        divisor = self.delta**l if self.is_boot else self.p**l
        q = self.q // divisor
        return self.__class__(
            self.b.divide(divisor, q), self.a.divide(divisor, q), self.is_boot
        )

    # Ciphertext bootstrapping status management

    def boot_to_nonboot(self):
        """
        Transform a ciphertext with is_boot = True to is_boot = False.

        Returns:
            CKKS:
                Transformed ciphertext with is_boot = False.
        """
        if not self.is_boot:
            return self

        l = floor(self.l * log(self.delta, self.p) + 1)
        q = self.moduli[l]
        return self.__class__(self.b % q, self.a % q, is_boot=False)

    def nonboot_to_boot(self):
        """
        Transform a ciphertext with is_boot = False to is_boot = True, all
        while also lifting it to the highest bootstrapping level.

        Returns:
            CKKS:
                Transformed ciphertext with is_boot = True.
        """
        if self.is_boot:
            return self

        q = self.moduli_boot[-1]
        return self.__class__(self.b.lift(q), self.a.lift(q), is_boot=True)

    def _check_boot_status(self, other):
        """
        Check if two ciphertexts have the same bootstrapping status.

        Args:
            other (CKKS):
                The other ciphertext to compare with.

        Raises:
            ValueError:
                If the ciphertexts do not have the same is_boot status.
        """
        if self.is_boot != other.is_boot:
            raise ValueError(
                "The ciphertexts do not have the same bootstrapping status."
            )

    # Addition and subtraction

    def __add__(self, other):
        """
        Add another ciphertext or plaintext.

        Args:
            other (CKKS or Poly or int):
                The ciphertext or plaintext to add.

        Returns:
            CKKS:
                Resulting ciphertext.
        """
        if isinstance(other, (int, Integer, Poly)):
            return self.__class__(self.b + other, self.a, self.is_boot)
        self._check_boot_status(other)
        return self.__class__(self.b + other.b, self.a + other.a, self.is_boot)

    def __radd__(self, other):
        """
        Add to another ciphertext or plaintext.

        Args:
            other (CKKS or Poly or int):
                The ciphertext or plaintext to add to.

        Returns:
            CKKS:
                Resulting ciphertext.
        """
        return self + other

    def __iadd__(self, other):
        """
        Add another ciphertext or plaintext, in-place.

        Args:
            other (CKKS or Poly or int):
                The other ciphertext or plaintext to add.

        Returns:
            CKKS:
                Resulting ciphertext.
        """
        return self + other

    def __neg__(self):
        """
        Negate the ciphertext.

        Returns:
            CKKS:
                Negated ciphertext.
        """
        return self.__class__(-self.b, -self.a, self.is_boot)

    def __sub__(self, other):
        """
        Subtract another ciphertext or plaintext.

        Args:
            other (CKKS or Poly or int):
                The ciphertext or plaintext to subtract.

        Returns:
            CKKS:
                Resulting ciphertext.
        """
        if isinstance(other, (int, Integer, Poly)):
            return self.__class__(self.b - other, self.a, self.is_boot)
        self._check_boot_status(other)
        return self.__class__(self.b - other.b, self.a - other.a, self.is_boot)

    def __rsub__(self, other):
        """
        Subtract from another ciphertext or plaintext.

        Args:
            other (CKKS or Poly or int):
                The ciphertext or plaintext to subtract from.

        Returns:
            CKKS:
                Resulting ciphertext.
        """
        if isinstance(other, (int, Integer)):
            return -(self - other)
        return other - self

    def __isub__(self, other):
        """
        Subtract another ciphertext or plaintext, in-place.

        Args:
            other (CKKS or Poly or int):
                The ciphertext or plaintext to subtract.

        Returns:
            CKKS:
                Resulting ciphertext.
        """
        return self - other

    # Multiplication

    def __mul__(self, other):
        """
        Multiply with another ciphertext or plaintext.

        Args:
            other (CKKS or Poly or int or 1j):
                The ciphertext or plaintext to multiply with.

        Returns:
            CKKS:
                Resulting ciphertext.
        """
        self._check_key_gen()
        if isinstance(other, (int, Integer, Poly)):
            return self.__class__(other * self.b, other * self.a, self.is_boot)
        if other == 1j:
            monomial = Poly.get_monomial(self.N // 2, self.N, self.q)
            return self * monomial

        self._check_boot_status(other)
        d0 = self.b * other.b
        d1 = self.a * other.b + other.a * self.b
        d2 = self.a * other.a
        d2_lift = d2.lift(d2.q * self.P)
        return self.__class__(
            d0 + (d2_lift * self.evk[0]).divide(self.P, d2.q),
            d1 + (d2_lift * self.evk[1]).divide(self.P, d2.q),
            self.is_boot,
        )

    def __rmul__(self, other):
        """
        Right multiply with another ciphertext or plaintext.

        Args:
            other (CKKS or Poly or int):
                The ciphertext or plaintext to right-multiply with.

        Returns:
            CKKS:
                Resulting ciphertext.
        """
        return self * other

    def __imul__(self, other):
        """
        Multiply with another ciphertext or plaintext, in-place.

        Args:
            other (CKKS or Poly or int):
                The ciphertext or plaintext to multiply with.

        Returns:
            CKKS:
                Resulting ciphertext.
        """
        return self * other

    def __matmul__(self, other):
        """
        Multiply and rescale with another ciphertext or plaintext.

        Args:
            other (CKKS or Poly or int):
                The ciphertext or plaintext to multiply with.

        Returns:
            CKKS:
                Resulting ciphertext.
        """
        return (self * other).rescale()

    def __rmatmul__(self, other):
        """
        Right multiply and rescale with another ciphertext or plaintext.

        Args:
            other (CKKS or Poly or int):
                The ciphertext or plaintext to right-multiply with.

        Returns:
            CKKS:
                Resulting ciphertext.
        """
        return self @ other

    def __imatmul__(self, other):
        """
        Multiply and rescale with another ciphertext or plaintext, in-place.

        Args:
            other (CKKS or Poly or int):
                The ciphertext or plaintext to multiply with.

        Returns:
            CKKS:
                Resulting ciphertext.
        """
        return self @ other

    # Applying polynomials

    def get_powers(self, d):
        """
        Compute all powers of the ciphertext self up to a given degree d.

        Args:
            d (int):
                Maximum exponent to compute.

        Returns:
            list:
                List of ciphertexts, where the i-th element is self raised to
                the power i = 0, ..., d. Powers are with respect to the
                operation @. For instance, squaring is self @ self.
        """
        powers = [0] * (d + 1)
        c = floor(log(d, 2))

        powers[0] = self.enc_poly_without_error(self.p)
        powers[1] = self

        for i in range(0, c):
            pow_two = 2**i
            powers[2 * pow_two] = powers[pow_two] @ powers[pow_two]
            for j in range(1, min(2 * pow_two, d + 1 - 2 * pow_two)):
                powers[2 * pow_two + j] = powers[2 * pow_two] @ powers[j]

        return powers

    def apply_poly(self, poly_coeffs):
        """
        Apply a polynomial to the ciphertext by minimizing the multiplicative
        depth.

        Args:
            poly_coeffs (list):
                List of coefficients (Poly instances or integers). These
                coefficients are assumed to be multiplied by the scaling factor
                p or delta already.

        Returns:
            CKKS:
                Resulting ciphertext.
        """
        if len(poly_coeffs) == 0:
            return self.enc_poly_without_error(0, self.is_boot)
        ct_out = self.enc_poly_without_error(poly_coeffs[0], self.is_boot)
        if len(poly_coeffs) == 1:
            return ct_out
        d = len(poly_coeffs) - 1
        powers = self.get_powers(d)
        for i in range(1, d + 1):
            ct_out += poly_coeffs[i] @ powers[i]
        return ct_out

    # Galois

    def key_switch(self, swk, P=None):
        """
        Switch keys using the provided switching key.

        Args:
            swk (tuple):
                Switching key (swk0, swk1).
            P (int, optional):
                Extra factor added to the modulus. Defaults to cls.P.

        Returns:
            CKKS:
                Resulting ciphertext.
        """
        self._check_key_gen()
        if P is None:
            P = self.P
        b = self.b
        a = self.a
        a_lift = a.lift(a.q * P)
        return self.__class__(
            b + (a_lift * swk[0]).divide(P, a.q),
            (a_lift * swk[1]).divide(P, a.q),
            self.is_boot,
        )

    def galois(self, k, swk=None):
        """
        Apply a Galois automorphism to the ciphertext.

        Args:
            k (int):
                Exponent of the Galois automorphism.
            swk (tuple, optional):
                Switching key for the automorphism.

        Returns:
            CKKS:
                Resulting ciphertext.

        Raises:
            KeyError:
                If the required switching key has not been generated.
        """
        k = k % (2 * self.N)
        if k == 1:
            return self
        if swk is None:
            if k not in self.galois_swk_dict:
                raise KeyError(
                    "Corresponding switching key swk not yet generated by "
                    "Alice."
                )
            swk = self.galois_swk_dict[k]
        return self.__class__(
            self.b.galois(k), self.a.galois(k), self.is_boot
        ).key_switch(swk)

    def conjugate(self, swk=None):
        """
        Apply complex conjugation to the ciphertext.

        Args:
            swk (tuple, optional):
                Switching key for conjugation.

        Returns:
            CKKS:
                Resulting ciphertext.
        """
        return self.galois(-1, swk)

    def rotate(self, k, n=None, swk=None):
        """
        Rotate the underlying plaintext vector by k slots.

        Args:
            k (int):
                Number of slots to rotate.
            n (int, optional):
                Length of the underlying plaintext vector.
            swk (tuple, optional):
                Switching key for the rotation.

        Returns:
            CKKS:
                Resulting ciphertext after rotation.
        """
        if n is None:
            n = self.n
        k = k % n
        return self.galois(5**k, swk)

    def conj_rotate(self, k, n=None, swk=None):
        """
        Conjugate and rotate the underlying plaintext vector by k slots.

        Args:
            k (int):
                Number of slots to rotate.
            n (int, optional):
                Length of the underlying plaintext vector.
            swk (tuple, optional):
                Switching key for the conjugation with rotation.

        Returns:
            CKKS:
                Resulting ciphertext after conjugation and rotation.
        """
        if n is None:
            n = self.n
        k = k % n
        return self.galois(-(5**k), swk)

    # Trace and product operations

    def trace(self, a, b, swks=None):
        """
        Perform a trace operation on the underlying plaintext vector, from a
        slots to b slots. Slot i of the resulting vector will be the sum of the
        slots i, i + b, i + 2 * b, ... of the original vector.

        Args:
            a (int):
                Length of the original underlying plaintext vector.
            b (int):
                Length of the resulting vector.
            swks (list, optional):
                List of switching keys for the required rotations.

        Returns:
            CKKS:
                Resulting ciphertext after applying the partial trace opertor.

        Raises:
            ValueError:
                If a is not a power of two.
            ValueError:
                If b does not divide a.
        """
        if a & (a - 1) != 0:
            raise ValueError(
                "The length a of the original vector is a power of two."
            )
        if a % b != 0:
            raise ValueError(
                "The length b of the resulting vector should divide "
                "the length a of the original vector."
            )

        log_a_over_b = int(round(log(a // b, 2)))
        if swks is None:
            swks = [None] * log_a_over_b
        for l in range(log_a_over_b):
            self += self.rotate(a // 2 ** (l + 1), n=a, swk=swks[l])
        return self

    def product(self, a, b, swks=None):
        """
        Perform a product operation on the underlying plaintext vector, from a
        slots to b slots. Slot i of the resulting vector will be the product of
        the slots i, i + b, i + 2 * b, ... of the original vector.

        Args:
            a (int):
                Length of the original underlying plaintext vector.
            b (int):
                Length of the resulting plaintext vector.
            swks (list, optional):
                List of switching keys for the required rotations.

        Returns:
            CKKS:
                Resulting ciphertext after applying the partial product
                opertor.

        Raises:
            ValueError:
                If a is not a power of two.
            ValueError:
                If b does not divide a.
        """
        if a & (a - 1) != 0:
            raise ValueError(
                "The length a of the original vector is a power of two."
            )
        if a % b != 0:
            raise ValueError(
                "The length b of the resulting vector should divide "
                "the length a of the original vector."
            )
        log_a_over_b = int(round(log(a // b, 2)))
        if swks is None:
            swks = [None] * log_a_over_b
        for l in range(log_a_over_b):
            self @= self.rotate(a // 2 ** (l + 1), swks[l])
        return self

    # Multiplication by special matrices

    @classmethod
    def get_poly_matrix(cls, A, is_boot=False):
        """
        Encode diagonals of matrix A as polynomials.

        Args:
            A (Multidiags):
                Square matrix whose diagonals are to be encoded.
            is_boot (bool, optional):
                Whether the encoded matrix is used during bootstrapping.
                Defaults to False.

        Returns:
            dict:
                Dictionary mapping diagonal indices to encoded polynomials.

        Raises:
            ValueError:
                If size of A does not divide N / 2.
        """
        if cls.N // 2 % A.n != 0:
            raise ValueError("The matrix size should divide N / 2.")

        return {
            i: cls.encode(A.get_diag(i), is_boot)
            for i in A.get_symm_diag_indices()
        }

    def BSGS_left_mult(self, poly_matrix):
        """
        Perform fast matrix-ciphertext multiplication using the Baby-Step
        Giant-Step algorithm.

        Args:
            poly_matrix (dict):
                Dictionary mapping diagonal indices to polynomials.

        Returns:
            CKKS:
                Resulting ciphertext.

        Raises:
            ValueError:
                If the matrix is not regular.
        """
        U = list(poly_matrix.keys())
        if 0 not in U:
            raise ValueError("The matrix is not regular.")
        if len(U) == 1:
            return self @ poly_matrix[0]
        d = U[1] - U[0]
        u_min = min(U)
        u_max = max(U)
        t = (u_max - u_min) // d + 1

        while t & t - 1 != 0:
            # Ensure that t is a power of two
            t += 1
        half_log_t = log(t, 2) / 2
        k0 = 2 ** floor(half_log_t)
        k1 = 2 ** ceil(half_log_t)

        rotated_ct = {}

        ct_out = self.enc_poly_without_error(0, self.is_boot)
        for i in range(k0):
            ct = self.enc_poly_without_error(0, self.is_boot)
            is_non_zero_term = False
            a = d * i * k1
            for j in range(k1):
                b = d * j + u_min
                c = a + b
                if c > u_max:
                    # Zero diagonal
                    continue
                if c in poly_matrix:
                    is_non_zero_term = True
                    if b not in rotated_ct:
                        rotated_ct[b] = self.rotate(b)
                    my_poly = poly_matrix[c]
                    my_poly = my_poly.galois(5 ** ((-a) % self.n))
                    ct += my_poly * rotated_ct[b]
            if is_non_zero_term:
                ct = ct.rotate(a)
                ct_out += ct

        return ct_out.rescale()

    @classmethod
    def get_BSGS_rotation_indices(cls, poly_matrix):
        """
        Calculate the rotation indices required for the fast matrix ciphertext
        multiplication using the Baby-Step Giant-Step algorithm. This is useful
        to precompute the corresponding switching keys.

        Args:
            poly_matrix (dict):
                Dictionary mapping diagonal indices to polynomials.

        Returns:
            list:
                The list of rotation indices.
        """
        U = list(poly_matrix.keys())
        if 0 not in U:
            raise ValueError("The matrix is not regular.")
        if len(U) == 1:
            return [0]
        d = U[1] - U[0]
        u_min = min(U)
        u_max = max(U)
        t = (u_max - u_min) // d + 1

        while t & t - 1 != 0:
            # Ensure that t is a power of two
            t += 1
        half_log_t = log(t, 2) / 2
        k0 = 2 ** floor(half_log_t)
        k1 = 2 ** ceil(half_log_t)

        rotation_indices = [(d * j + u_min) % cls.n for j in range(k1)] + [
            (d * i * k1) % cls.n for i in range(k0)
        ]

        return rotation_indices

    #  Bootstrapping

    @classmethod
    def config_bootstrap(cls, sk, d=7, r=8, log_radix=1, print_messages=False):
        """
        Perform precomputations needed for bootstrapping.

        Args:
            sk (Poly):
                Secret key.
            d (int, optional):
                Degree of the Taylor polynomial used for EvalMod. Defaults to
                7.
            r (int, optional):
                Number of successive squarings at the end of EvalMod. Defaults
                to 8.
            log_radix (int, optional):
                Grouping parameter controlling the trade-off between
                multiplicative depth and number of multiplications during
                bootstrapping. Defaults to 1; this is fast, but consumes many
                CKKS levels.
            print_messages (bool, optional):
                Whether to print information messages. Defaults to False.

        """
        log_radix = min(log_radix, log(cls.n, 2))

        cls._check_key_gen()

        if print_messages:
            print("Preparing EvalMod coefficients...")

        cls.r = r
        scalar = cls.delta * np.pi * 1j / (2**r * cls.q0)
        cls.EvalMod_encoded_coeffs = [
            cls.encode(scalar**k / factorial(k), True) for k in range(d + 1)
        ]

        if print_messages:
            print(
                "Generating and encoding matrices required for CoeffToSlot "
                "and SlotToCoeff..."
            )

        if cls.n != 1:
            grouped_F = get_grouped_F(
                cls.n, log_radix, False
            ).copy()  # Matrices for StC
            grouped_iF = get_grouped_F(
                cls.n, log_radix, True
            ).copy()  # Matrices for CtS

        else:
            grouped_F = [Multidiags.get_identity(1)]
            grouped_iF = [Multidiags.get_identity(1)]

        # Division by N / (2 * n) incorporated in the CtS matrices
        k = len(grouped_iF)
        scalar = np.power(2 * cls.n / cls.N, 1 / k)
        for i in range(k):
            grouped_iF[i] = grouped_iF[i].incorporate_scalar(scalar)

        # Product by q0 / (4 * delta * pi) incorporated in the StC matrices
        k = len(grouped_F)
        scalar = np.power(cls.q0 / (4 * cls.delta * np.pi), 1 / k)
        for i in range(k):
            grouped_F[i] = grouped_F[i].incorporate_scalar(scalar)

        grouped_poly_F = [
            cls.get_poly_matrix(A, is_boot=True) for A in grouped_F
        ]
        grouped_poly_iF = [
            cls.get_poly_matrix(A, is_boot=True) for A in grouped_iF
        ]
        cls.grouped_poly_F = (grouped_poly_F, grouped_poly_iF)

        if print_messages:
            print("Generating missing switching keys...")

        rotation_indices = []
        for poly_matrix in grouped_poly_F + grouped_poly_iF:
            rotation_indices += cls.get_BSGS_rotation_indices(poly_matrix)
        for k in rotation_indices:
            cls.get_galois_swk(5**k, sk)
        cls.get_galois_swk(-1, sk)  # For conjugation
        for k in range(log(cls.N, 2)):
            # For rotation by powers of two
            cls.get_galois_swk(5 ** (2**k), sk)

        if print_messages:
            print("The bootstrapping configuration is done!")

    @classmethod
    def _check_config_bootstrap(cls):
        """
        Check if the scheme is configured for bootstrapping.

        Raises:
            RuntimeError:
                 If config_bootstrap needs to be called first.
        """
        try:
            cls.r
        except AttributeError:
            raise RuntimeError(
                "Ask Alice to configure the bootstrapping first."
            )

    def bootstrap(self):
        """
        Refresh the ciphertext self by increasing its level.

        Returns:
            CKKS:
                A new ciphertext with increased level.
        """

        self._check_config_bootstrap()

        # ModRaise

        ct = (self % self.q0).nonboot_to_boot()
        ct = ct.trace(self.N // 2, self.n)

        # CoeffToSlot

        grouped_poly_iF = self.grouped_poly_F[1]
        for A in grouped_poly_iF:
            ct = ct.BSGS_left_mult(A)
        ct_conj = ct.conjugate()
        ct0 = ct + ct_conj
        ct1 = 1j * (ct_conj - ct)
        cts = [ct0, ct1]

        # EvalMod

        for i in range(2):
            ct = cts[i].apply_poly(self.EvalMod_encoded_coeffs)
            for _ in range(self.r):
                ct = ct @ ct
            cts[i] = 1j * (ct.conjugate() - ct)

        # SlotToCoeff

        ct = cts[0] + 1j * cts[1]
        grouped_poly_F = self.grouped_poly_F[0]
        for A in grouped_poly_F:
            ct = ct.BSGS_left_mult(A)

        return ct.boot_to_nonboot()

    def get_precision(self, other, sk, n=None):
        """
        Compute the average number of bits of precision preserved in the
        coefficients of the underlying plaintext polynomial of other, relative
        to self.

        Args:
            other (CKKS):
                The ciphertext to compare against self.
            sk (Poly):
                The secret key.
            n (int, optional):
                Number of slots to consider. Defaults to self.n.

        Returns:
            float:
                The average number of bits of precision preserved from self to
                other.
        """
        if n is None:
            n = self.n

        pt = self.dec_to_poly(sk).lift(0)
        e = pt - other.dec_to_poly(sk).lift(0)

        coeffs_pt = pt.get_symmetric_coeffs()
        coeffs_pt = np.array(
            [abs(int(coeffs_pt[i * self.N // (2 * n)])) for i in range(2 * n)]
        )

        coeffs_e = e.get_symmetric_coeffs()
        coeffs_e = np.array(
            [abs(int(coeffs_e[i * self.N // (2 * n)])) for i in range(2 * n)]
        )

        bits_pt = np.log2(np.abs(coeffs_pt) + 1)
        bits_e = np.log2(np.abs(coeffs_e) + 1)

        precision = np.round(np.mean(bits_pt - bits_e), 2)
        if precision < 0:
            return 0
        return precision

    # Representation

    def __repr__(self):
        """
        Return a string representation of the ciphertext. Provides basic
        information including degree, modulus, and level.

        Returns:
            str:
                String representation of the ciphertext.
        """
        if self.is_boot:
            str = [
                f"A CKKS ciphertext for bootstrapping ",
                f"with degree N = 2^{log(self.N, 2)} ",
                f"and modulus q = ",
                f"(2^{log(self.p*self.q0,2)}) * ",
                f"(2^{log(self.delta,2)})^{self.l}",
                f" (level {self.l} out of {self.L_boot})",
            ]
        else:
            str = [
                f"A CKKS ciphertext ",
                f"with degree N = 2^{log(self.N, 2)} ",
                f"and modulus q = ",
                f"(2^{log(self.q0,2)}) * (2^{log(self.p,2)})^{self.l}",
                f" (level {self.l} out of {self.L})",
            ]
        return "".join(str)
