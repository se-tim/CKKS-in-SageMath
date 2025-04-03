from sage.all import *
from .lattice_estimator.estimator import *
from .poly import *
from .fast_dft import *


class CKKS:
    """
    CKKS Homomorphic Encryption Scheme.

    This class implements the CKKS homomorphic encryption scheme, including key
    generation, encryption, decryption, homomorphic operations, and supports
    bootstrapping.
    """

    @classmethod
    def config(cls, N, n, L_boot, q0, p, delta, print_messages=True):
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
                Smallest modulus.
            p (int):
                Scaling factor outside of bootstrapping.
            delta (int):
                Scaling factor during bootstrapping.

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
        cls.q0 = q0  # Smallest modulus
        cls.p = p  # Scaling factor outside of bootstrapping
        cls.delta = delta  # Scaling factor during bootstrapping

        cls.L = floor(
            L_boot * log(delta, p) + 1
        )  # Maximal level outside of bootstrapping
        cls.L_boot = L_boot  # Maximal level for bootstrapping

        cls.moduli = [
            cls.q0 * cls.p**i for i in range(cls.L + 1)
        ]  # The ct moduli used outside of bootstrapping
        cls.moduli_boot = [
            cls.p * cls.q0 * cls.delta**i for i in range(cls.L_boot + 1)
        ]  # The ct moduli used during bootstrapping

        # Dictionary containing the polynomial versions of the grouped
        # matrices F_{2 * n, l} and iF_{2 * n, l}
        cls.grouped_poly_F_dict = {}

        cls.CtS_StC_poly_dict = (
            {}
        )  # Dictionary containing polynomials required for CtS and StC

        if print_messages:
            print("The CKKS configuration is done!\n")

    # Key generation

    @classmethod
    def key_gen(cls, h=None, sk=None, P=None, sigma=None, print_messages=True):
        """
        Generate the secret key, public key and evaluation key for the scheme.

        Args:
            h (int, optional):
                Hamming weight of the secret key. Defaults to
                2 ** (log(N, 2) // 2 - 1) if neither h nor sk are passed.
            sk (Poly, optional):
                Secret key for scheme. Defaults to a ternary polynomial of
                Hamming weight equal to h.
            P (int, optional):
                Factor for the evaluation key modulus. Defaults to the biggest
                modulus accepted by the scheme.
            sigma (float, optional):
                Standard deviation for error polynomials. Defaults to 3.2.

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
            cls.h = 2 ** (log(cls.N, 2) // 2 - 1)

        cls.P = cls.moduli_boot[-1] if P is None else P
        cls.sigma = 3.2 if sigma is None else sigma

        q_evk = cls.P * cls.moduli_boot[-1]
        if sk is not None:
            cls.sk = sk % q_evk
        else:
            cls.sk = Poly.get_random_ternary1(cls.N, q_evk, cls.h)

        e = Poly.get_random_normal(cls.N, cls.moduli_boot[-1], cls.sigma)
        pk1 = Poly.get_random_uniform(cls.N, cls.moduli_boot[-1])
        pk0 = e - pk1 * cls.sk
        cls.pk = (pk0, pk1)

        e = Poly.get_random_normal(cls.N, q_evk, cls.sigma)
        evk1 = Poly.get_random_uniform(cls.N, q_evk)
        evk0 = e + cls.P * cls.sk**2 - evk1 * cls.sk
        cls.evk = (evk0, evk1)

        cls.galois_swk_dict = {}  # Dictionary containing the swk for Galois

        # Check if the secret key is ternary; else it is binary
        if sk is None:
            sk_is_ternary = True
        else:
            sk_is_ternary = False
            for i in range(cls.N):
                if cls.sk.coeffs[i] == -1:
                    sk_is_ternary = True

        if print_messages:
            print("The key generation is done!")
            try:
                security = cls.get_security(sk_is_binary=not (sk_is_ternary))
                print(
                    f"Estimated security: 2^({log(security,2).n(digits=3)}) "
                    f"operations."
                )
            except (RuntimeError, TypeError):
                print("Estimated security: very low.")
            print()

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
    def get_swk(cls, sk_x, sk):
        """
        Generate a switching key to switch from secret key sk_x to sk.

        Args:
            sk_x (Poly):
                Secret key to switch from.
            sk (Poly):
                Secret key to switch to.

        Returns:
            tuple:
                Switching key (swk0, swk1).
        """
        cls._check_key_gen()
        q_evk = cls.P * cls.moduli_boot[-1]
        e = Poly.get_random_normal(cls.N, q_evk, cls.sigma)
        swk1 = Poly.get_random_uniform(cls.N, q_evk)  # Slow
        swk0 = e + cls.P * sk_x - swk1 * sk
        return (swk0, swk1)

    @classmethod
    def get_galois_swk(cls, k, sk=None):
        """
        Get switching key for applying Galois automorphism.

        Args:
            k (int):
                Exponent of the Galois automorphism.
            sk (Poly, optional):
                Secret key. Required if the switching key is not already
                generated.

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
        swk = cls.get_swk(sk_x, sk)
        cls.galois_swk_dict[k] = swk
        return swk

    # Security estimation

    @classmethod
    def get_security(cls, sk_is_binary=False):
        """
        Use the LWE estimator to compute the security level based on the
        current CKKS parameters.

        Args:
            sk_is_binary (bool):
                Whether the secret key is binary. Defaults to False.

        Returns:
            float:
                Logarithm (base 2) of estimated number of operations required
                to break the scheme.
        """
        (my_plus, my_minus) = (
            (cls.h // 2, cls.h // 2) if sk_is_binary == False else (cls.h, 0)
        )  # Number of coefficients in sk equal to +1 and -1, respectively

        params = LWE.Parameters(
            cls.N,
            cls.moduli_boot[-1] * cls.P,
            Xs=ND.SparseTernary(n=cls.N, p=my_plus, m=my_minus),
            Xe=ND.DiscreteGaussian(n=cls.N, stddev=cls.sigma),
        )
        return LWE.primal_usvp(params)["rop"]

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
        Add another ciphertext or plaintext, in place.

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
        Subtract another ciphertext or plaintext, in place.

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
                The ciphertext or plaintext to right multiply with.

        Returns:
            CKKS:
                Resulting ciphertext.
        """
        return self * other

    def __imul__(self, other):
        """
        Multiply with another ciphertext or plaintext, in place.

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
                The ciphertext or plaintext to right multiply with.

        Returns:
            CKKS:
                Resulting ciphertext.
        """
        return self @ other

    def __imatmul__(self, other):
        """
        Multiply and rescale with another ciphertext or plaintext, in place.

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

    def key_switch(self, swk):
        """
        Switch keys using the provided switching key.

        Args:
            swk (tuple):
                Switching key (swk0, swk1).

        Returns:
            CKKS:
                Resulting ciphertext.
        """
        self._check_key_gen()
        b = self.b
        a = self.a
        a_lift = a.lift(a.q * self.P)
        return self.__class__(
            b + (a_lift * swk[0]).divide(self.P, a.q),
            (a_lift * swk[1]).divide(self.P, a.q),
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

    # Trace and product operators

    def trace(self, a, b, swks=None):
        """
        Perform a partial trace operation on the underlying plaintext vector,
        from a slots to b slots. Slot i of the resulting vector will be the
        sum of the slots i, i + b, i + 2 * b, ... of the original vector.

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

        log_a_over_b = log(a / b, 2)
        if swks is None:
            swks = [None] * log_a_over_b
        for l in range(log_a_over_b):
            self += self.rotate(a // 2 ** (l + 1), n=a, swk=swks[l])
        return self

    def product(self, a, b, swks=None):
        """
        Perform a partial product operation on the underlying plaintext vector,
        from a slots to b slots. Slot i of the resulting vector will be the
        product of the slots i, i + b, i + 2 * b, ... of the original vector.

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

        log_a_over_b = log(a / b, 2)
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
        Perform fast matrix ciphertext multiplication using the Baby Step Giant
        Step algorithm.

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
        multiplication using the Baby Step Giant Step algorithm. This is useful
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
    def config_bootstrap(cls, sk, log_radix=1, print_messages=True):
        """
        Perform precomputations needed for bootstrapping.

        Args:
            sk (Poly):
                Secret key.
            log_radix (int, optional):
                Grouping parameter controlling the trade off between
                multiplicative depth and number of multiplications during
                bootstrapping. Defaults to 1; this is fast, but consumes many
                CKKS levels.
        """
        if 2**log_radix > cls.n:
            log_radix = log(cls.n, 2)

        cls._check_key_gen()

        if not cls.CtS_StC_poly_dict and 2 * cls.n != cls.N:
            if print_messages:
                print("Encoding relevant vectors as polynomials...")
            cls.CtS_StC_poly_dict["one_and_i"] = cls.encode(
                np.array(
                    [(1 - 1j) * (i < cls.n) + 1j for i in range(2 * cls.n)]
                ),
                is_boot=True,
            )
            cls.CtS_StC_poly_dict["one_and_zero"] = cls.encode(
                np.array([0.5 * (i < cls.n) for i in range(2 * cls.n)]),
                is_boot=True,
            )
            cls.CtS_StC_poly_dict["i_and_zero"] = cls.encode(
                np.array([-0.5 * 1j * (i >= cls.n) for i in range(2 * cls.n)]),
                is_boot=True,
            )

        if log_radix not in cls.grouped_poly_F_dict:
            if print_messages:
                print(
                    "Generating matrices required for CoeffToSlot and "
                    "SlotToCoeff..."
                )
            grouped_F = get_grouped_F(cls.n, log_radix, False)
            grouped_iF = get_grouped_F(cls.n, log_radix, True)

            if print_messages:
                print("Encoding these matrices as polynomials...")
            grouped_poly_F = [
                cls.get_poly_matrix(A, is_boot=True) for A in grouped_F
            ]
            grouped_poly_iF = [
                cls.get_poly_matrix(A, is_boot=True) for A in grouped_iF
            ]
            cls.grouped_poly_F_dict[log_radix] = (
                grouped_poly_F,
                grouped_poly_iF,
            )
        else:
            grouped_poly_F, grouped_poly_iF = cls.grouped_poly_F_dict[
                log_radix
            ]

        if print_messages:
            print("Generating missing switching keys...")

        rotation_indices = [] if 2 * cls.n == cls.N else [cls.n]
        for poly_matrix in grouped_poly_F + grouped_poly_iF:
            rotation_indices += cls.get_BSGS_rotation_indices(poly_matrix)

        for k in rotation_indices:
            # For CtS and StC
            cls.get_galois_swk(5**k, sk)

        cls.get_galois_swk(-1, sk)  # For conjugation

        for k in range(log(cls.N, 2)):
            # For trace operator
            cls.get_galois_swk(5 ** (2**k), sk)

        if print_messages:
            print("The bootstrapping configuration is done!\n")

    @classmethod
    def _check_config_bootstrap(cls, log_radix):
        """
        Check if the scheme is configured for bootstrapping with the specified
        grouping parameter s.

        Args:
            log_radix (int):
                Grouping parameter controlling the trade off between
                multiplicative depth and number of multiplications during
                bootstrapping. Defaults to 1; this is fast, but consumes many
                CKKS levels.

        Raises:
            RuntimeError:
                 If config_bootstrap needs to be called first.
        """
        if 2**log_radix > cls.n:
            log_radix = log(cls.n, 2)
        if log_radix not in cls.grouped_poly_F_dict:
            raise RuntimeError(
                "Ask Alice to configure the bootstrapping first."
            )

    def CoeffToSlot(self, log_radix=1):
        """
        Perform homomorphic encoding (coefficients to slots) with respect to
        the scaling factor delta, to obtain ciphertext(s) whose plaintext
        vector(s) encrypt(s) the coefficients of the initial plaintext
        polynomial.

        Args:
            log_radix (int, optional):
                Grouping parameter controlling the trade off between
                multiplicative depth and number of multiplications during
                CoeffToSlot. Defaults to 1; this is fast, but consumes many
                CKKS levels.

        Returns:
            list:
                - If we have full slots (n = N / 2), returns a two element list
                  [ct0, ct1] of ciphertexts.
                - Otherwise, returns a list [ct] with a single ciphertext.

        """
        if 2**log_radix > self.n:
            log_radix = log(self.n, 2)

        self._check_config_bootstrap(log_radix)

        grouped_poly_iF = self.grouped_poly_F_dict[log_radix][1]

        ct = self.nonboot_to_boot()
        for A in grouped_poly_iF:
            ct = ct.BSGS_left_mult(A)

        ct_conj = ct.conjugate()

        if 2 * self.n == self.N:
            # Full slots
            ct0 = self.delta // 2 * (ct + ct_conj)
            ct1 = -self.delta // 2 * (1j * (ct - ct_conj))
            return [ct0, ct1]

        ct0 = ct + ct_conj
        ct1 = ct - ct_conj
        return [
            ct0 * self.CtS_StC_poly_dict["one_and_zero"]
            + ct1 * self.CtS_StC_poly_dict["i_and_zero"]
        ]

    @classmethod
    def SlotToCoeff(cls, cts, log_radix=1):
        """
        Perform homomorphic decoding (slots to coefficients) with respect to
        the scaling factor delta, to obtain a ciphertext whose underlying
        plaintext polynomial has as coefficients the values contained in the
        plaintext vector(s) of the input ciphertext(s), multiplied by the
        factor delta.

        Args:
            cts (tuple):
                - If we have full slots (n = N / 2), cts is a list [ct0, ct1]
                  of ciphertexts.
                - Otherwise, cts is a list [ct] of a single ciphertext.
            log_radix (int, optional):
                Grouping parameter controlling the trade off between
                multiplicative depth and number of multiplications during
                SlotToCoeff. Defaults to 1; this is fast, but consumes many
                CKKS levels.

        Returns:
            CKKS:
                The desired ciphertext.
        """
        if 2**log_radix > cls.n:
            log_radix = log(cls.n, 2)

        cls._check_config_bootstrap(log_radix)

        grouped_poly_F = cls.grouped_poly_F_dict[log_radix][0]

        if 2 * cls.n == cls.N:
            # Full slots
            ct0, ct1 = cts
            ct = ct0 + 1j * ct1

        else:
            ct = cts[0]
            ct @= cls.CtS_StC_poly_dict["one_and_i"]
            ct += ct.rotate(cls.n, 2 * cls.n)

        for A in grouped_poly_F:
            ct = ct.BSGS_left_mult(A)

        return ct

    def EvalMod(self, M, d, r):
        """
        Reduce the ciphertext modulo a positive integer M by using a degree d
        polynomial to approximate the complex exponential function, which is
        then successively squared r times to improve precision, before the
        imaginary part is taken.

        Args:
            M (int):
                The modulus to reduce the ciphertext by.
            d (int):
                Degree of the polynomial used for approximation.
            r (int):
                Number of times the polynomial approximation is squared.

        Returns:
            CKKS:
                A new ciphertext approximating self reduced modulo M.
        """
        scaling_factor = self.delta if self.is_boot else self.p
        f0 = self.encode(
            2 * scaling_factor * np.pi * 1j / (2**r * M), self.is_boot
        )
        ct = (f0 @ self).rescale()  # Multiply by 2 * pi * 1j / (2**r * M)
        encoded_coeffs = [
            self.encode(1 / factorial(k), self.is_boot) for k in range(d + 1)
        ]
        ct = ct.apply_poly(encoded_coeffs)
        for _ in range(r):
            ct = ct @ ct
        ct = ct - ct.conjugate()
        f1 = self.encode(M / (4 * np.pi * 1j))
        return f1 @ ct

    def bootstrap(self, log_radix=1, d=7, r=7):
        """
        Refresh the ciphertext self by increasing its level.

        Args:
            log_radix (int, optional):
                Grouping parameter controlling the trade off between
                multiplicative depth and number of multiplications during
                bootstrapping. Defaults to 1; this is fast, but consumes many
                CKKS levels.
            d (int, optional):
                Degree of the polynomial used in EvalMod. Defaults to 7.
            r (int, optional):
                Number of times the polynomial approximation is squared in
                EvalMod. Defaults to 7.

        Returns:
            CKKS:
                A new ciphertext with increased level.
        """
        if 2**log_radix > self.n:
            log_radix = log(self.n, 2)

        self._check_config_bootstrap(log_radix)

        ct = (self % self.q0).nonboot_to_boot()

        if 2 * self.n < self.N:
            ct = ct.trace(self.N // 2, self.n)
            ct = (
                (self.delta * 2 * self.n // self.N) * ct
            ).rescale()  # Divide by N / (2 * n)

        cts = ct.CoeffToSlot(log_radix)  # List of one or two ciphertexts

        for i in range(len(cts)):
            cts[i] = cts[i].EvalMod(self.q0, d, r)

        ct = self.SlotToCoeff(cts, log_radix)
        ct = (self.delta // self.p) @ ct  # Divide by p
        return ct.boot_to_nonboot()

    def get_poly_precision(self, other, sk, n=None):
        """
        Compute the average number of bits of precision preserved in the
        coefficients of the underlying plaintext polynomial of other, relative
        to self.

        Args:
            other (CKKS):
                The ciphertext to compare against self.
            sk (Poly):
                Secret key of scheme.
            n (int, optional):
                Number of slots to consider. Defaults to self.n.

        Returns:
            float:
                The average number of bits of precision preserved from self to
                other.
        """
        if n is None:
            n = self.n

        pt = self.dec_to_poly(sk)
        e = pt - other.dec_to_poly(sk)

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

        return np.mean(bits_pt - bits_e)

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
                f" (level {self.l} out of {self.L_boot}).",
            ]
        else:
            str = [
                f"A CKKS ciphertext ",
                f"with degree N = 2^{log(self.N, 2)} ",
                f"and modulus q = ",
                f"(2^{log(self.q0,2)}) * (2^{log(self.p,2)})^{self.l}",
                f" (level {self.l} out of {self.L}).",
            ]
        return "".join(str)
