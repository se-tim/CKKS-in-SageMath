from estimator import *

load("poly.py")
load("fast_dft.py")


class CKKS:
    """
    CKKS Homomorphic Encryption Scheme.

    This class implements the CKKS homomorphic encryption scheme, including key
    generation, encryption, decryption, homomorphic operations, and supports
    bootstrapping.
    """

    @classmethod
    def config(cls, N, n, L, q0, delta):
        """
        Configure basic CKKS parameters.

        Args:
            N (int):
                Ring degree (must be a power of two).
            n (int):
                Number of slots (must be a strict divisor of N).
            L (int):
                Maximal level.
            q0 (int):
                Smallest modulus.
            delta (int):
                Scaling factor.

        Raises:
            ValueError:
                If N is not a power of two.
            ValueError:
                If n is not a strict divisor of N.
        """
        print("Setting the scheme's parameters...")
        if N & N - 1 != 0:
            raise ValueError("N must be a power of two.")
        if N % n != 0 or n >= N:
            raise ValueError("n must be a strict divisor of N.")
        cls.N = N  # Ring degree
        cls.n = n  # Number of slots
        cls.L = L  # Maximal level
        cls.q0 = q0  # Smallest modulus
        cls.delta = delta  # Scaling factor
        cls.moduli = [
            cls.q0 * cls.delta**i for i in range(L + 1)
        ]  # All ct moduli for the scheme

        print("Generating matrices required for encoding and decoding...")
        get_grouped_E(n, 1, inverse=True)
        get_grouped_E(n, 1, inverse=False)

        # Dictionary containing the polynomial versions of the groups of
        # matrices F_{2n, l} and iF_{2n, l}
        cls.grouped_poly_F_dict = {}

        cls.CtS_StC_poly_dict = (
            {}
        )  # Dictionary containing polynomials required for CtS and StC

        print("The CKKS configuration is done!\n")

    # Key generation

    @classmethod
    def key_gen(cls, h=None, P=None, sigma=3.2):
        """
        Generate the secret key, public key and evaluation key for the scheme.

        Args:
            h (int, optional):
                Hamming weight of the secret key. Defaults to
                2 ** (log(N, 2) // 2 - 1).
            P (int, optional):
                Factor for the evaluation key modulus. Defaults to the biggest
                modulus accepted by the scheme, namely cls.moduli[-1].
            sigma (float, optional):
                Standard deviation for error polynomials. Defaults to 3.2.

        Raises:
            RuntimeError:
                If CKKS is not configured before key generation.
        """
        print("Setting the scheme's key parameters...")
        try:
            cls.N
        except AttributeError:
            raise RuntimeError("Ask Alice to configure CKKS first.")
        cls.h = (
            2 ** (log(cls.N, 2) // 2 - 1) if h is None else h
        )  # Hamming weight of secret key
        cls.P = (
            cls.moduli[-1] if P is None else P
        )  # Required for modulus of evaluation key
        cls.sigma = sigma  # Standard deviation for error polynomials

        print("Generating secret key...")
        q_evk = cls.P * cls.moduli[-1]
        cls.sk = Poly.get_random_ternary1(cls.N, q_evk, cls.h)  # Secret key

        print("Generating public key...")
        e = Poly.get_random_normal(cls.N, cls.moduli[-1], cls.sigma)
        pk1 = Poly.get_random_uniform(cls.N, cls.moduli[-1])
        pk0 = e - pk1 * cls.sk
        cls.pk = (pk0, pk1)  # Public key

        print("Generating evaluation key...")
        e = Poly.get_random_normal(cls.N, q_evk, cls.sigma)
        evk1 = Poly.get_random_uniform(cls.N, q_evk)
        evk0 = e + cls.P * cls.sk**2 - evk1 * cls.sk
        cls.evk = (evk0, evk1)  # Evaluation key

        cls.galois_swk_dict = {}  # Dictionary containing the swk for Galois

        try:
            print(
                f"Estimated security: "
                f"2^({log(cls.get_security(),2).n(digits=3)})."
            )
        except (RuntimeError, TypeError):
            print("Estimated security: very low.")
        print("The key generation is done!\n")

    @classmethod
    def _check_key_gen(cls):
        """
        Check if CKKS is configured and keys have been generated.

        Raises:
            RuntimeError:
                If CKKS is not configured or keys have not been generated.
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
        q_evk = cls.P * cls.moduli[-1]
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
            TypeError:
                If sk is None and the switching key is not already generated.
        """
        k = k % (2 * cls.N)
        cls._check_key_gen()
        if k in cls.galois_swk_dict:
            return cls.galois_swk_dict[k]
        if sk is None:
            raise TypeError("Secret key sk required.")
        sk_x = sk.galois(k)
        swk = cls.get_swk(sk_x, sk)
        cls.galois_swk_dict[k] = swk
        return swk

    # Security estimation

    @classmethod
    def get_security(cls):
        """
        Use the LWE estimator to compute the security level based on the
        current CKKS parameters.

        Returns:
            float:
                Estimated security level.
        """
        params = LWE.Parameters(
            cls.N,
            cls.moduli[-1] * cls.P,
            Xs=ND.SparseTernary(n=cls.N, p=cls.h // 2, m=cls.h // 2),
            Xe=ND.DiscreteGaussian(n=cls.N, stddev=cls.sigma),
        )
        return LWE.primal_usvp(params)["rop"]

    # Initialization of ciphertexts

    @classmethod
    def _check_modulus(cls, q):
        """
        Check if the current parameters support the given modulus q.

        Args:
            q (int): Modulus to check.

        Raises:
            ValueError:
                If modulus q is not supported by the current CKKS parameters.
        """
        if q != 0 and q not in cls.moduli:
            raise ValueError(
                "Modulus q not supported by current CKKS parameters."
            )

    def __init__(self, b, a):
        """
        Initialize a CKKS ciphertext with tuple (b, a).

        Args:
            b (Poly):
                Polynomial b.
            a (Poly):
                Polynomial a.

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
        if a.q == 0:
            q = self.moduli[-1]
            self.b, self.a = b % q, a % q
            self.q = q  # Modulus of ciphertext
            self.l = self.L  # Level of ciphertext
        else:
            self._check_modulus(a.q)
            self.b, self.a = b, a
            self.q = a.q
            self.l = self.moduli.index(a.q)  # Level of ciphertext

    # Encoding and decoding

    @classmethod
    def encode(cls, z):
        """
        Encode a complex vector z (in bit reversed order) into an integer
        polynomial. This process requires applying the inverse of a variant of
        the DFT matrix, and scaling the result by the factor delta to achieve
        the necessary precision.

        Args:
            z (np.ndarray or np.complex128):
                A complex vector (of length 1 if z is an np.complex128) to
                encode. The vector is expected to be in bit reversed order.

        Returns:
            Poly:
                Integer polynomial encoding the vector z, but in bit reversed
                order.

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

        if n == 1:
            # Single value, representing a vector whose entries are
            # identical to this value.
            a = int(round(cls.delta * real(z)))
            b = int(round(cls.delta * imag(z)))
            return Poly.get_constant(
                a, cls.N, cls.moduli[-1]
            ) + b * Poly.get_monomial(cls.N // 2, cls.N, cls.moduli[-1])

        grouped_iE = get_grouped_E(n, 1, inverse=True)

        w = copy(z)
        for A in grouped_iE:
            w = A.BSGS_mult(w)
        w *= cls.delta

        pt0_temp = np.round(real(w)).astype(int)
        pt1_temp = np.round(imag(w)).astype(int)
        pt0, pt1 = [np.zeros(cls.N // 2).astype(int) for _ in range(2)]

        l = cls.N // (2 * n)
        for i in range(n):
            pt0[l * i] = pt0_temp[i]
            pt1[l * i] = pt1_temp[i]

        return Poly(np.concatenate((pt0, pt1)), cls.N, cls.moduli[-1])

    @classmethod
    def decode(cls, pt, n=None):
        """
        Decode an integer polynomial to a complex vector (in bit reversed
        order) of length n. This process is effectively the inverse of the
        previously described encode method.

        Args:
            pt (Poly):
                Integer polynomial to decode.
            n (int, optional):
                Length of the complex vector. Defaults to cls.n.

        Returns:
            np.ndarray:
                Complex vector z of length n (in bit reversed order).

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
        w /= cls.delta

        return w

    # Encryption and decryption

    @classmethod
    def enc_poly_with_sk(cls, pt, sk):
        """
        Encrypt a plaintext polynomial with the secret key.

        Args:
            pt (Poly):
                Plaintext polynomial to encrypt.
            sk (Poly):
                Secret key.

        Returns:
            CKKS:
                Encrypted ciphertext.

        Raises:
            ValueError:
                If the modulus of pt is not supported.
        """
        cls._check_modulus(pt.q)
        q = cls.moduli[-1] if pt.q == 0 else pt.q
        e = Poly.get_random_normal(cls.N, q, cls.sigma)
        a = Poly.get_random_uniform(cls.N, q)
        b = pt + e - a * sk
        return cls(b, a)

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

        Raises:
            ValueError:
                If the modulus of pt is not supported.
        """
        cls._check_modulus(pt.q)
        q = cls.moduli[-1] if pt.q == 0 else pt.q
        v = Poly.get_random_ternary2(cls.N, q)
        e0, e1 = [
            Poly.get_random_normal(cls.N, q, cls.sigma) for _ in range(2)
        ]
        return cls(pt + v * pk[0] + e0, v * pk[1] + e1)

    @classmethod
    def enc_poly_without_error(cls, pt):
        """
        Encrypt a plaintext polynomial without error (insecure).

        Args:
            pt (Poly or int):
                Plaintext polynomial or integer to encrypt.

        Returns:
            CKKS:
                Encrypted ciphertext.
        """
        if isinstance(pt, (int, Integer)):
            a = Poly.get_constant(0, cls.N, cls.moduli[-1])
            b = Poly.get_constant(pt, cls.N, cls.moduli[-1])
            return cls(b, a)
        cls._check_modulus(pt.q)
        q = cls.moduli[-1] if pt.q == 0 else pt.q
        a = Poly.get_constant(0, cls.N, q)
        b = pt
        return cls(b, a)

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

    # Error computation

    def get_noise(self, z, sk):
        """
        Calculate the exact noise in the ciphertext.

        Args:
            z (np.ndarray):
                Original complex vector that was encrypted.
            sk (Poly):
                Secret key.

        Returns:
            float:
                Error between the original vector and the decrypted vector
                (infinity norm of difference).
        """
        pt = self.dec_to_poly(sk)
        w = CKKS.decode(pt)
        return np.max(np.abs(z - w))

    def get_precision_bits(self, z, sk):
        """
        Calculate the number of correct bits in the ciphertext compared to
        the original plaintext vector.

        Args:
            z (np.ndarray):
                Original plaintext vector.
            sk (Poly):
                Secret key.

        Returns:
            float:
                The total number of correct bits in the real and imaginary
                parts combined.
        """
        w = z - CKKS.decode(self.dec_to_poly(sk))
        num_real_bits = -np.log(
            np.max(np.abs(np.real(w))) / np.max(np.abs(np.real(z)))
        ) / np.log(2)
        num_imag_bits = -np.log(
            np.max(np.abs(np.imag(w))) / np.max(np.abs(np.imag(z)))
        ) / np.log(2)
        return num_real_bits + num_imag_bits

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
        return self.__class__(self.b % q, self.a % q)

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
        if q is None:
            q = self.moduli[-1]
        return __class__(self.b.lift(q), self.a.lift(q))

    def rescale(self, l=1):
        """
        Rescale by dividing by delta**l.

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
        divisor = self.delta**l
        q = self.q // divisor
        return self.__class__(
            self.b.divide(divisor, q), self.a.divide(divisor, q)
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
            return self.__class__(self.b + other, self.a)
        return self.__class__(self.b + other.b, self.a + other.a)

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
        return self.__class__(-self.b, -self.a)

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
            return self.__class__(self.b - other, self.a)
        return self.__class__(self.b - other.b, self.a - other.a)

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

        Raises:
            RuntimeError:
                If keys have not been generated.
        """
        self._check_key_gen()
        if isinstance(other, (int, Integer, Poly)):
            return self.__class__(other * self.b, other * self.a)
        if self == 1j:
            monomial = Poly.get_monomial(other.N // 2, other.N, other.q)
            return other * monomial
        if other == 1j:
            monomial = Poly.get_monomial(self.N // 2, self.N, self.q)
            return self * monomial
        d0 = self.b * other.b
        d1 = self.a * other.b + other.a * self.b
        d2 = self.a * other.a
        d2_lift = d2.lift(d2.q * self.P)
        return self.__class__(
            d0 + (d2_lift * self.evk[0]).divide(self.P, d2.q),
            d1 + (d2_lift * self.evk[1]).divide(self.P, d2.q),
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

        powers[0] = self.enc_poly_without_error(self.delta)
        powers[1] = self

        for i in range(0, c):
            p = 2**i
            powers[2 * p] = powers[p] @ powers[p]
            for j in range(1, min(2 * p, d + 1 - 2 * p)):
                powers[2 * p + j] = powers[2 * p] @ powers[j]

        return powers

    def apply_poly(self, poly_coeffs):
        """
        Apply a polynomial to the ciphertext by minimizing the multiplicative
        depth.

        Args:
            poly_coeffs (list):
                List of coefficients (Poly instances or integers). These
                coefficients are assumed to be multiplied by the scaling factor
                delta already.

        Returns:
            CKKS:
                Resulting ciphertext.
        """
        if len(poly_coeffs) == 0:
            return self.enc_poly_without_error(0)
        ct_out = self.enc_poly_without_error(poly_coeffs[0])
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
        return self.__class__(self.b.galois(k), self.a.galois(k)).key_switch(
            swk
        )

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
        Rotate the underlying plaintext vector (in bit reversed order) by k
        slots.

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

    # Multiplication by special matrices

    @classmethod
    def get_poly_matrix(cls, A):
        """
        Apply bit reversal to the diagonals of matrix A, then encode them as
        polynomials.

        Args:
            A (Multidiags):
                Square matrix of size n x n whose diagonals are to be encoded.

        Returns:
            dict:
                Dictionary mapping diagonal indices to encoded polynomials.

        Raises:
            ValueError:
                If matrix A has incorrect size (must be n x n).
        """

        def bit_rev_vector(z, num_bits):
            # Put coefficients of vector into bit reversed order
            return np.array([z[bit_rev(i, num_bits)] for i in range(len(z))])

        if A.n != cls.n:
            raise ValueError("Matrix is of wrong size.")

        poly_matrix = {}
        indices = A.get_symm_diag_indices()

        num_bits = log(cls.n, 2)
        for i in indices:
            z = bit_rev_vector(A.get_diag(i), num_bits)
            poly_matrix[i] = cls.encode(z)

        return poly_matrix

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
        t = len(U)
        if t == 1:
            return self @ poly_matrix[0]
        d = U[1] - U[0]
        u_min = min(U)
        u_max = max(U)

        while t & t - 1 != 0:
            # Ensure that t is a power of two
            t += 1
        half_log_2_t = log(t, 2) / 2
        k0 = 2 ** floor(half_log_2_t)
        k1 = 2 ** ceil(half_log_2_t)

        rotated_ct = {
            d * j + u_min: self.rotate(d * j + u_min) for j in range(k1)
        }

        ct_out = self.enc_poly_without_error(0)
        for i in range(k0):
            ct0 = self.enc_poly_without_error(0)
            a = d * i * k1
            for j in range(k1):
                b = d * j + u_min
                c = a + b
                if c > u_max:
                    # Zero diagonal
                    continue
                p = poly_matrix[c]
                p = p.galois(5 ** ((-a) % self.n))
                ct1 = p @ rotated_ct[b]
                ct0 = ct0 + ct1
            ct0 = ct0.rotate(a)
            ct_out += ct0
        return ct_out

    #  Bootstrapping

    @classmethod
    def config_bootstrap(cls, sk, s=1):
        """
        Perform precomputations needed for bootstrapping.

        Args:
            sk (Poly):
                Secret key.
            s (int, optional):
                Grouping parameter controlling the trade off between
                multiplicative depth and number of multiplications during
                bootstrapping. Defaults to 1; this is fast, but consumes many
                CKKS levels.
        """
        if 2**s > cls.n:
            s = log(cls.n, 2)

        cls._check_key_gen()

        if not cls.CtS_StC_poly_dict and 2 * cls.n != cls.N:
            print("Encoding relevant vectors as polynomials...")
            cls.CtS_StC_poly_dict["one_and_i"] = CKKS.encode(
                np.array(
                    [(1 - 1j) * (i % 2 == 0) + 1j for i in range(2 * cls.n)]
                )
            )
            cls.CtS_StC_poly_dict["one_and_zero"] = CKKS.encode(
                np.array([0.5 * (i % 2 == 0) for i in range(2 * cls.n)])
            )
            cls.CtS_StC_poly_dict["i_and_zero"] = CKKS.encode(
                np.array([-0.5 * 1j * (i % 2 == 1) for i in range(2 * cls.n)])
            )

        if s not in cls.grouped_poly_F_dict:
            print(
                "Generating matrices required for CoeffToSlot and "
                "SlotToCoeff..."
            )
            grouped_F = get_grouped_F(cls.n, s, False)
            grouped_iF = get_grouped_F(cls.n, s, True)

            print("Encoding these matrices as polynomials...")
            grouped_poly_F = [cls.get_poly_matrix(A) for A in grouped_F]
            grouped_poly_iF = [cls.get_poly_matrix(A) for A in grouped_iF]

            cls.grouped_poly_F_dict[s] = (
                grouped_poly_F,
                grouped_poly_iF,
            )
        else:
            grouped_poly_F, grouped_poly_iF = cls.grouped_poly_F_dict[s]

        print("Generating missing switching keys...")

        rotation_indices = [] if 2 * cls.n == cls.N else [cls.n]
        for poly_matrix in grouped_poly_F + grouped_poly_iF:
            U = list(poly_matrix.keys())
            t = len(U)
            d = U[1] - U[0]
            u_min = min(U)

            while t & t - 1 != 0:
                # Ensure that t is a power of two
                t += 1
            half_log_2_t = log(t, 2) / 2
            k0 = 2 ** floor(half_log_2_t)
            k1 = 2 ** ceil(half_log_2_t)

            rotation_indices += [
                (d * j + u_min) % cls.n for j in range(k1)
            ] + [(d * i * k1) % cls.n for i in range(k0)]

        for k in rotation_indices:
            # For CtS and StC
            cls.get_galois_swk(5**k, sk)

        cls.get_galois_swk(-1, sk)  # For conjugation

        l = log(cls.N // (2 * cls.n), 2)
        for k in range(l):
            # For partial sum
            cls.get_galois_swk(5 ** (cls.n * 2**k), sk)

        print("The bootstrapping configuration is done!\n")

    @classmethod
    def _check_config_bootstrap(cls, s):
        """
        Check if the scheme is configured for bootstrapping with the specified
        grouping parameter s.

        Args:
            s (int):
                Grouping parameter controlling the trade off between
                multiplicative depth and number of multiplications during
                bootstrapping. Defaults to 1; this is fast, but consumes many
                CKKS levels.

        Raises:
            RuntimeError:
                 If config_bootstrap needs to be called first.
        """
        if 2**s > cls.n:
            s = log(cls.n, 2)
        if s not in cls.grouped_poly_F_dict:
            raise RuntimeError(
                "Ask Alice to configure the bootstrapping first."
            )

    def CoeffToSlot(self, s=1):
        """
        Perform homomorphic encoding (coefficients to slots) to obtain
        ciphertext(s) whose plaintext vector(s) encrypt(s) the coefficients of
        the initial plaintext polynomial.

        Args:
            s (int, optional):
                Grouping parameter controlling the trade off between
                multiplicative depth and number of multiplications during
                CoeffToSlot. Defaults to 1; this is fast, but consumes many
                CKKS levels.

        Returns:
            tuple or CKKS:
                - If we have full slots (n = N / 2), returns a two element list
                  [ct0, ct1] of ciphertexts.
                - Otherwise, returns a list [ct] with a single ciphertext.

        Raises:
            RuntimeError:
                If config_bootstrap needs to be called first.
        """
        if 2**s > self.n:
            s = log(self.n, 2)

        self._check_config_bootstrap(s)

        grouped_poly_iF = self.grouped_poly_F_dict[s][1]

        ct = copy(self)
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
    def SlotToCoeff(cls, cts, s=1):
        """
        Perform homomorphic decoding (slots to coefficients) to obtain a
        ciphertext whose underlying plaintext polynomial has as coefficients
        the values contained in the plaintext vector(s) of the input
        ciphertext(s), multiplied by the factor delta.

        Args:
            cts (tuple):
                - If we have full slots (n = N / 2), cts is a list [ct0, ct1]
                  of ciphertexts.
                - Otherwise, cts is a list [ct] of a single ciphertext.
            s (int, optional):
                Grouping parameter controlling the trade off between
                multiplicative depth and number of multiplications during
                SlotToCoeff. Defaults to 1; this is fast, but consumes many
                CKKS levels.

        Returns:
            CKKS:
                The desired ciphertext.

        Raises:
            RuntimeError:
                If config_bootstrap needs to be called first.
        """
        if 2**s > cls.n:
            s = log(cls.n, 2)

        cls._check_config_bootstrap(s)

        grouped_poly_F = cls.grouped_poly_F_dict[s][0]

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

    def partial_sum(self):
        """
        Return a ciphertext corresponding to summing all slots of the plaintext
        vector whose indices differ by a multiple of n.

        Returns:
            CKKS:
                A new ciphertext where each slot contains the sum of specific
                slots from the original plaintext vector.
        """
        ct = copy(self)
        l = log(self.N // (2 * self.n), 2)
        for k in range(l):
            ct_x = ct.rotate(self.n * 2**k, self.N // 2)
            ct += ct_x
        return ct

    def EvalMod(self, p, d, r):
        """
        Reduce the ciphertext modulo a positive integer p by using a degree d
        polynomial to approximate the complex exponential function, which is
        then successively squared r times to improve precision, before the
        imaginary part is taken.

        Args:
            p (int):
                The modulus to reduce the ciphertext by.
            d (int):
                Degree of the polynomial used for approximation.
            r (int):
                Number of times the polynomial approximation is squared.

        Returns:
            CKKS:
                A new ciphertext approximating self reduced modulo p.
        """
        f0 = self.encode(2 * self.delta * np.pi * 1j / (2**r * p))
        ct0 = (f0 @ self).rescale()
        encoded_coeffs = [self.encode(1 / factorial(k)) for k in range(d + 1)]
        ct1 = ct0.apply_poly(encoded_coeffs)
        for _ in range(r):
            ct1 = ct1 @ ct1
        ct1 = ct1 - ct1.conjugate()
        f1 = self.encode(p / (4 * np.pi * 1j))
        return f1 @ ct1

    def bootstrap(self, s=1, d=7, r=5):
        """
        Refresh the ciphertext self by increasing its level.

        Args:
            s (int, optional):
                Grouping parameter controlling the trade off between
                multiplicative depth and number of multiplications during
                bootstrapping. Defaults to 1; this is fast, but consumes many
                CKKS levels.
            d (int, optional):
                Degree of the polynomial used in EvalMod. Defaults to 7.
            r (int, optional):
                Number of times the polynomial approximation is squared in
                EvalMod. Defaults to 5.

        Returns:
            CKKS:
                A new ciphertext with increased level.

        Raises:
            RuntimeError:
                If config_bootstrap needs to be called first.
        """
        if 2**s > self.n:
            s = log(self.n, 2)

        self._check_config_bootstrap(s)

        ct = (self % self.q0).lift()

        if self.n < self.N / 2:
            ct = ct.partial_sum()
            ct = (
                (self.delta * 2 * self.n // self.N) * ct
            ).rescale()  # Remove factor N // (2 * n)

        cts = ct.CoeffToSlot(s)  # List of one or two ciphertexts

        for i in range(len(cts)):
            cts[i] = (
                (self.delta**2 // self.q0) * (cts[i].rescale())
            ).rescale()  # Divide by q0
            cts[i] = cts[i].EvalMod(1, d, r)
            cts[i] = (self.q0 // self.delta) * cts[
                i
            ]  # Multiply by q0 // delta

        ct = CKKS.SlotToCoeff(cts, s)  # Also removes factor 1 / delta
        return ct

    # Representation

    def __repr__(self):
        """
        Return a string representation of the ciphertext. Provides basic
        information including degree, modulus, and level.

        Returns:
            str:
                String representation of the ciphertext.
        """
        str = [
            f"A CKKS ciphertext with degree N = 2^{log(self.N, 2)} ",
            f"and modulus q = ",
            f"(2^{log(self.q0,2)}) * (2^{log(self.delta,2)})^{self.l}",
            f" (level {self.l} out of {self.L})",
        ]
        return "".join(str)
