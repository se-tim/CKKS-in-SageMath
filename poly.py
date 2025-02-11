import numpy as np
from sage.libs.ntl import *  # Fast polynomial multiplication
from sage.libs.ntl.ntl_ZZ_p import (
    ntl_ZZ_p_random_element,
)  # Random elements from ntl_ZZ_p


class Poly:
    """
    Class representing integer polynomials modulo X ** N + 1 and an integer
    modulus q.
    """

    # Some static class methods

    galois_int_dic = (
        {}
    )  # Dictionary containing the permutation of integers (see below)

    @classmethod
    def _get_galois_integers(cls, k, N):
        """
        Multiply the integers in [0, N - 1] by k and consider them modulo
        2 * N.

        Args:
            k (int):
                Exponent for the Galois automorphism.
            N (int):
                Degree of the polynomial ring.

        Returns:
            list:
                A list of tuples (j, b), where j is the index and b is a bit
                indicating sign.

        Raises:
            ValueError:
                If k is not coprime to 2 * N.
        """
        if (k, N) in cls.galois_int_dic:
            return cls.galois_int_dic[(k, N)]
        if gcd(k, 2 * N) != 1:
            raise ValueError(
                "The exponent k must be coprime to 2 and to the degree N."
            )
        indices = []
        for i in range(N):
            j = (k * i) % (2 * N)
            if j >= N:
                indices.append((j - N, 1))
            else:
                indices.append((j, 0))
        cls.galois_int_dic[(k, N)] = indices
        return indices

    @staticmethod
    def _get_ntl_coeffs(coeffs, N, q):
        """
        Create NTL polynomial coefficients from coeffs.

        Args:
            coeffs (list, np.ndarray or NTL coefficients):
                Coefficients.
            N (int):
                Degree of the polynomial ring.
            q (int):
                Modulus for coefficients.

        Returns:
            NTL polynomial coefficients.
        """
        if isinstance(coeffs, list):
            return ntl.ZZX(coeffs) if q == 0 else ntl.ZZ_pX(coeffs, q)
        if isinstance(coeffs, np.ndarray):
            coeffs = list(coeffs)
            return ntl.ZZX(coeffs) if q == 0 else ntl.ZZ_pX(coeffs, q)
        input_q = (
            0
            if isinstance(coeffs, ntl_ZZX.ntl_ZZX)
            else coeffs.get_modulus_context().modulus()
        )
        if input_q == q:
            return coeffs
        if q == 0:
            return ntl.ZZX([coeffs[i].lift_centered() for i in range(N)])
        if input_q == 0:
            return ntl.ZZ_pX(coeffs.list(), q)
        return coeffs.convert_to_modulus(ntl.ZZ_pContext(q))

    @classmethod
    def _get_constant_ntl_coeffs(cls, k, N, q):
        """
        Create NTL coefficients which are all equal to k.

        Args:
            k (int):
                Constant value for coefficients.
            N (int):
                Degree of the polynomial ring.
            q (int):
                Modulus for coefficients.

        Returns:
            NTL polynomial coefficients with all coefficients equal to k.
        """
        return Poly._get_ntl_coeffs([k] * N, N, q)

    @staticmethod
    def _quotient(coeffs, N):
        """
        Compute coeffs modulo X ** N + 1.

        Args:
            coeffs:
                NTL polynomial coefficients.
            N (int):
                Degree of the polynomial ring.

        Returns:
            NTL polynomial coefficients modulo X ** N + 1.
        """
        if coeffs.degree() >= N:
            return coeffs.truncate(N) - coeffs.right_shift(N)
        return coeffs

    # Initialization

    def __init__(self, coeffs, N, q=0):
        """
        Initialize a polynomial modulo q and modulo X ** N + 1 from coeffs.

        Args:
            coeffs (list, np.ndarray or NTL coefficients):
                Coefficients.
            N (int):
                Degree of the polynomial ring.
            q (int, optional):
                Modulus for coefficients. Defaults to 0.
        """
        self.coeffs = self._quotient(self._get_ntl_coeffs(coeffs, N, q), N)
        self.q = q
        self.N = N

    # Methods for internal use

    def _check_degree(self, N):
        """
        Raise an error if self.N is different from the input N.

        Args:
            N (int):
                Degree to check against.

        Raises:
            ValueError:
                If degrees are incompatible.
        """
        if self.N != N:
            raise ValueError("Incompatible degrees N.")

    def _check_modulus(self, q):
        """
        Raise an error if q does not divide self.q.

        Args:
            q (int):
                Modulus to check.

        Raises:
            ValueError:
                If moduli are incompatible.
        """
        if self.q % q != 0:
            raise ValueError("Incompatible moduli q.")

    # Symmetric coefficients

    def get_symmetric_coeffs(self):
        """
        Return self.coeffs as a list, each coefficient between
        (-self.q / 2, self.q / 2], unless q == 0.

        Returns:
            list:
                List of symmetric coefficients.
        """
        if self.q == 0:
            return [self.coeffs[i] for i in range(self.N)]
        return [self.coeffs[i].lift_centered() for i in range(self.N)]

    # Modular reduction and lifting

    def __mod__(self, q):
        """
        Compute self % q.

        Args:
            q (int):
                Modulus.

        Returns:
            Poly:
                Resulting polynomial modulo q.
        """
        self._check_modulus(q)
        if self.q == q:
            return self
        return self.__class__(
            self._get_ntl_coeffs(self.coeffs, self.N, q), self.N, q
        )

    def __imod__(self, q):
        """
        Compute self %= q.

        Args:
            q (int):
                Modulus.

        Returns:
            Poly:
                Resulting polynomial modulo q.
        """
        return self % q

    def lift(self, q):
        """
        Lift the polynomial to a multiple of the original modulus.

        Args:
            q (int):
                New modulus (must be a multiple of self.q).

        Returns:
            Poly:
                Lifted polynomial.

        Raises:
            ValueError:
                If q is not a multiple of self.q.
        """
        if q % self.q != 0:
            raise ValueError(
                "Can only lift to a multiple of the original modulus."
            )
        q_half = self.q // 2
        if q == 0:
            return self.__class__(self.get_symmetric_coeffs(), self.N, q)
        r = self._get_constant_ntl_coeffs(q_half, self.N, self.q)
        s = self._get_constant_ntl_coeffs(q_half, self.N, q)
        coeffs = self.coeffs.left_shift(0)
        coeffs += r
        coeffs = ntl.ZZ_pX([coeffs[i] for i in range(self.N)], q) - s
        return self.__class__(coeffs, self.N, q)

    # Addition and subtraction

    def _add_sub(self, other, sign):
        """
        Internal method for addition and subtraction.

        Args:
            other (Poly or int):
                The other polynomial or integer.
            sign (int):
                0 for addition, 1 for subtraction.

        Returns:
            Poly:
                Result of the addition or subtraction.
        """
        if isinstance(other, (int, Integer)):
            return self.__class__(
                self.coeffs
                + self._get_ntl_coeffs([(-1) ** sign * other], self.N, self.q),
                self.N,
                self.q,
            )
        if isinstance(other, Poly):
            self._check_degree(other.N)
            if self.q == other.q:
                return self.__class__(
                    (
                        self.coeffs + other.coeffs
                        if sign == 0
                        else self.coeffs - other.coeffs
                    ),
                    self.N,
                    self.q,
                )
            if self.q == 0:
                return self.__class__(
                    (
                        self._get_ntl_coeffs(self.coeffs, other.N, other.q)
                        + other.coeffs
                        if sign == 0
                        else self._get_ntl_coeffs(
                            self.coeffs, other.N, other.q
                        )
                        - other.coeffs
                    ),
                    other.N,
                    other.q,
                )
            if other.q == 0:
                return self.__class__(
                    (
                        self.coeffs
                        + self._get_ntl_coeffs(other.coeffs, self.N, self.q)
                        if sign == 0
                        else self.coeffs
                        - self._get_ntl_coeffs(other.coeffs, self.N, self.q)
                    ),
                    self.N,
                    self.q,
                )
            if self.q < other.q:
                other._check_modulus(self.q)
                return self.__class__(
                    (
                        self.coeffs
                        + self._get_ntl_coeffs(other.coeffs, self.N, self.q)
                        if sign == 0
                        else self.coeffs
                        - self._get_ntl_coeffs(other.coeffs, self.N, self.q)
                    ),
                    self.N,
                    self.q,
                )
            self._check_modulus(other.q)
            return self.__class__(
                (
                    self._get_ntl_coeffs(self.coeffs, other.N, other.q)
                    + other.coeffs
                    if sign == 0
                    else self._get_ntl_coeffs(self.coeffs, other.N, other.q)
                    - other.coeffs
                ),
                other.N,
                other.q,
            )
        if sign == 0:
            return other + self
        return -other + self

    def __add__(self, other):
        """
        Add another polynomial or integer.

        Args:
            other (Poly or int):
                The polynomial or integer to add.

        Returns:
            Poly:
                Resulting polynomial.
        """
        return self._add_sub(other, 0)

    def __radd__(self, other):
        """
        Add to another polynomial or integer.

        Args:
            other (Poly or int):
                The polynomial or integer to add to.

        Returns:
            Poly:
                Resulting polynomial.
        """
        return self + other

    def __iadd__(self, other):
        """
        Add another polynomial or integer, in place.

        Args:
            other (Poly or int):
                The polynomial or integer to add.

        Returns:
            Poly:
                Resulting polynomial.
        """
        return self + other

    def __neg__(self):
        """
        Negate the polynomial.

        Returns:
            Poly:
                Negated polynomial.
        """
        return self.__class__(-self.coeffs, self.N, self.q)

    def __sub__(self, other):
        """
        Subtract another polynomial or integer.

        Args:
            other (Poly or int):
                The polynomial or integer to subtract.

        Returns:
            Poly:
                Resulting polynomial.
        """
        return self._add_sub(other, 1)

    def __rsub__(self, other):
        """
        Subtract from another polynomial or integer.

        Args:
            other (Poly or int):
                The polynomial or integer to subtract from.

        Returns:
            Poly:
                Resulting polynomial.
        """
        if isinstance(other, (int, Integer)):
            return -(self._add_sub(other, 1))
        return other._add_sub(self, 1)

    def __isub__(self, other):
        """
        Subtract another polynomial or integer, in place.

        Args:
            other (Poly or int):
                The polynomial or integer to subtract.

        Returns:
            Poly:
                Resulting polynomial.
        """
        return self - other

    # Multiplication

    def __mul__(self, other):
        """
        Multiply with another polynomial or integer.

        Args:
            other (Poly or int):
                The polynomial or integer to multiply with.

        Returns:
            Poly:
                Resulting polynomial.
        """
        if isinstance(other, (int, Integer)):
            return self.__class__(
                self.coeffs * self._get_ntl_coeffs([other], self.N, self.q),
                self.N,
                self.q,
            )
        if isinstance(other, Poly):
            self._check_degree(other.N)
            if self.q == other.q:
                return self.__class__(
                    self.coeffs * other.coeffs,
                    self.N,
                    self.q,
                )
            if self.q == 0:
                return self.__class__(
                    self._get_ntl_coeffs(self.coeffs, other.N, other.q)
                    * other.coeffs,
                    other.N,
                    other.q,
                )
            if other.q == 0:
                return self.__class__(
                    self.coeffs
                    * self._get_ntl_coeffs(other.coeffs, self.N, self.q),
                    self.N,
                    self.q,
                )
            if self.q < other.q:
                other._check_modulus(self.q)
                return self.__class__(
                    self.coeffs
                    * self._get_ntl_coeffs(other.coeffs, self.N, self.q),
                    self.N,
                    self.q,
                )
            self._check_modulus(other.q)
            return self.__class__(
                self._get_ntl_coeffs(self.coeffs, other.N, other.q)
                * other.coeffs,
                other.N,
                other.q,
            )
        return other * self

    def __rmul__(self, other):
        """
        Right multiply with another polynomial or integer.

        Args:
            other (Poly or int):
                The polynomial or integer to right multiply with.

        Returns:
            Poly:
                Resulting polynomial.
        """
        return self * other

    def __imul__(self, other):
        """
        Multiply with another polynomial or integer, in place.

        Args:
            other (Poly or int):
                The polynomial or integer to multiply with.

        Returns:
            Poly:
                Resulting polynomial.
        """
        return self * other

    def __pow__(self, exp):
        """
        Take to the power of exp.

        Args:
            exp (int):
                Exponent (must be non negative).

        Returns:
            Poly:
                Resulting polynomial.

        Raises:
            ValueError:
                If exponent is negative.
        """
        if exp < 0:
            raise ValueError("Exponent must be non negative.")
        if exp == 0:
            return self.__class__([1], self.N, self.q)
        if exp == 1:
            return self
        if exp % 2 == 0:
            result = self ** (exp // 2)
            return result * result
        result = self ** ((exp - 1) // 2)
        return self * result * result

    def shift(self, k):
        """
        Multiply the polynomial by the polynomial X ** k.

        Args:
            k (int):
                Exponent k.

        Returns:
            Poly:
                Resulting polynomial.
        """
        k = k % (2 * self.N)
        coeffs = self.coeffs.left_shift(k)
        return self.__class__(coeffs, self.N, self.q)

    # Division

    def divide(self, divisor, q=None):
        """
        Divide the coefficients of the polynomial by divisor, then round.

        Args:
            divisor (int):
                Divisor (must be positive).
            q (int, optional):
                New modulus. Defaults to None.

        Returns:
            Poly:
                Resulting polynomial.

        Raises:
            ZeroDivisionError:
                If divisor is zero.
            ValueError:
                If divisor is negative.
        """
        if divisor <= 0:
            if divisor == 0:
                raise ZeroDivisionError("Divisor must be positive.")
            else:
                raise ValueError("Divisor must be positive.")
        if divisor == 1 and (q is None or q == self.q):
            return self
        if self.q == 0:
            coeffs = np.empty(self.N, dtype=object)
            for i in range(self.N):
                value = ZZ(self.coeffs[i])
                value //= divisor
                coeffs[i] = value
            return self.__class__(coeffs, self.N, self.q if q is None else q)
        elif q == self.q // divisor:  # Fastest
            divisor_half = divisor // 2
            coeffs = (
                self.coeffs
                + self._get_constant_ntl_coeffs(divisor_half, self.N, self.q)
            )._right_pshift(ntl.ZZ(divisor))
            return self.__class__(coeffs, self.N, q)
        else:
            coeffs = np.empty(self.N, dtype=object)
            q_half = self.q // 2
            for i in range(self.N):
                value = ZZ(self.coeffs[i])
                if value > q_half:
                    value -= self.q
                value //= divisor
                coeffs[i] = value
            return self.__class__(coeffs, self.N, self.q if q is None else q)

    # Galois

    def galois(self, k):
        """
        Apply a Galois automorphism: send X ** i to X ** {i * k}.

        Args:
            k (int):
                Exponent for the automorphism.

        Returns:
            Poly:
                Resulting polynomial.
        """
        two_N = 2 * self.N
        k = k % two_N
        if k == 1:
            return self
        if k == two_N - 1:
            coeffs = self.coeffs.left_shift(0)
            c0 = coeffs[0]
            coeffs[0] = 0
            coeffs[self.N] = 1
            coeffs = -coeffs.reverse()
            coeffs[0] = c0
            return self.__class__(coeffs, self.N, self.q)
        coeffs = self._get_ntl_coeffs([], self.N, self.q)
        permuted_indices = self._get_galois_integers(k, self.N)
        for i in range(self.N):
            j, b = permuted_indices[i]
            if b == 1:
                coeffs[j] = -self.coeffs[i]
            else:
                coeffs[j] = self.coeffs[i]
        return self.__class__(coeffs, self.N, self.q)

    # Special polynomials

    constant_poly_dic = {}  # Dictionary containing constant polynomials

    @classmethod
    def get_constant(cls, k, N, q):
        """
        Return a constant polynomial corresponding to the integer k.

        Args:
            k (int):
                Constant value.
            N (int):
                Degree of the polynomial ring.
            q (int):
                Modulus for coefficients.

        Returns:
            Poly:
                Constant polynomial.
        """
        if (k, N, q) in cls.constant_poly_dic:
            return cls.constant_poly_dic[(k, N, q)]
        a = cls([k], N, q)
        cls.constant_poly_dic[(k, N, q)] = a
        return a

    monomial_dic = {}  # Dictionary containing monomials

    @classmethod
    def get_monomial(cls, k, N, q):
        """
        Return the monomial corresponding to X ** k.

        Args:
            k (int):
                Exponent k.
            N (int):
                Degree of the polynomial ring.
            q (int):
                Modulus for coefficients.

        Returns:
            Poly:
                Monomial polynomial.
        """
        k = k % (2 * N)
        if (k, N, q) in cls.monomial_dic:
            return cls.monomial_dic[(k, N, q)]
        a = cls(
            cls._get_ntl_coeffs([1 if k < N else -1], N, q).left_shift(k % N),
            N,
            q,
        )
        cls.monomial_dic[(k, N, q)] = a
        return a

    @classmethod
    def get_random_uniform(cls, N, q):
        """
        Return a polynomial with coefficients chosen uniformly at random
        in [0, q - 1].

        Args:
            N (int):
                Degree of the polynomial ring.
            q (int):
                Modulus for coefficients.

        Returns:
            Poly:
                Random polynomial.

        Raises:
            ZeroDivisionError:
                If q is zero.
            ValueError:
                If q is negative.
        """
        if q == 0:
            raise ZeroDivisionError(
                "Cannot draw uniformly at random modulo 0."
            )
        if q < 0:
            raise ValueError("Modulus q must be a positive integer.")
        coeffs = [ntl_ZZ_p_random_element(q) for _ in range(N)]
        return cls(coeffs, N, q)

    @classmethod
    def get_random_ternary1(cls, N, q, h=None):
        """
        Return a polynomial with h coefficients chosen uniformly at random
        in [-1, 1].

        Args:
            N (int):
                Degree of the polynomial ring.
            q (int):
                Modulus for coefficients.
            h (int, optional):
                Number of non zero coefficients. Defaults to N - 1.

        Returns:
            Poly:
                Random ternary polynomial.

        Raises:
            ValueError:
                If h is not in [0, N - 1].
        """
        h = N - 1 if h is None else h
        if not (0 <= h < N):
            raise ValueError("h must lie in [0, N - 1].")
        k = np.random.randint(0, h + 1)
        coeffs = [1] * k + [-1] * (h - k) + [0] * (N - h)
        np.random.shuffle(coeffs)
        return cls(coeffs, N, q)

    @classmethod
    def get_random_ternary2(cls, N, q, rho=0.5):
        """
        Return a polynomial with coefficients from [-1, 1], where
        prob(0) = rho and prob(-1) = prob(1).

        Args:
            N (int):
                Degree of the polynomial ring.
            q (int):
                Modulus for coefficients.
            rho (float, optional):
                Probability of zero coefficient. Defaults to 0.5.

        Returns:
            Poly:
                Random ternary polynomial.

        Raises:
            ValueError:
                If rho is not in [0, 1].
        """
        if not (0 <= rho <= 1):
            raise ValueError("rho must lie in [0, 1].")
        coeffs = np.random.choice(
            range(-1, 2), size=N, p=[(1 - rho) / 2, rho, (1 - rho) / 2]
        )
        return cls(coeffs, N, q)

    @classmethod
    def get_random_normal(cls, N, q, sigma=3.2):
        """
        Return a polynomial with coefficients chosen normally around zero
        with standard deviation sigma.

        Args:
            N (int):
                Degree of the polynomial ring.
            q (int):
                Modulus for coefficients.
            sigma (float, optional):
                Standard deviation. Defaults to 3.2.

        Returns:
            Poly:
                Random polynomial with normal distribution.
        """
        coeffs = np.random.normal(0, sigma, N).round().astype(int)
        return cls(coeffs, N, q)

    # Representation

    def __repr__(self):
        """
        Represent the polynomial as a string.

        Returns:
            str:
                String representation of the polynomial.
        """
        terms = []
        coeffs = self.get_symmetric_coeffs()
        poly_str = ""
        for i in range(len(coeffs)):
            coeff = coeffs[i]
            if coeff != 0:
                if i > 0 and coeff == 1:
                    term = f"+ "
                elif i > 0 and coeff == -1:
                    term = f"- "
                else:
                    term = f"+ {coeff}" if coeff > 0 else f"- {-coeff}"
                if i == 1:
                    term += "X"
                elif i > 1:
                    term += f"X^{i}"
                terms.append(term)
        if not terms:
            terms.append("0")
        if self.q != 0:
            terms.append(f"mod({factor(self.q)})")
        poly_str = " ".join(terms) if terms else "0"
        if poly_str[:2] == "+ ":
            poly_str = poly_str[2:]
        return poly_str
