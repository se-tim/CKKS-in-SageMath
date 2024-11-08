import numpy as np
import pickle  # Saving in files

# Multidiags class


class Multidiags:
    """
    Class representing a matrix by its non zero diagonals.
    """

    def __init__(self, n, diags={}):
        """
        Initialize the Multidiags instance.

        Constructs an n x n matrix represented by its non zero diagonals
        provided as a dictionary.

        Args:
            n (int):
                The size of the matrix (n x n).
            diags (dict, optional):
                A dictionary where keys are diagonal indices and values are
                numpy arrays representing the diagonals. Defaults to {}.
        """
        self.n = n
        self.diags = {}
        for i in diags:
            self.diags[i % n] = diags[i]

    @classmethod
    def get_identity(cls, n):
        """
        Get the identity matrix of size n x n.

        Args:
            n (int):
                The size of the matrix.

        Returns:
            Multidiags:
                An identity matrix of size n x n.
        """
        return cls(n, {0: np.ones(n, dtype=np.complex128)})

    def get_diag_indices(self):
        """
        Get the indices of the non zero diagonals.

        Returns:
            list:
                List of diagonal indices.
        """
        return list(self.diags.keys())

    def get_symm_diag_indices(self):
        """
        Get the indices of the non zero diagonals, reduced symmetrically
        around 0.

        Returns:
            list:
                List of symmetric diagonal indices.
        """
        half_n = self.n // 2
        J = [j if j <= half_n else j - self.n for j in self.diags]
        J.sort()
        return J

    def get_diag(self, i):
        """
        Get the i-th diagonal of the matrix.

        Args:
            i (int):
                The index of the diagonal.

        Returns:
            np.ndarray:
                The i-th diagonal as a numpy array.
        """
        i %= self.n
        return (
            self.diags[i]
            if i in self.diags
            else np.zeros(self.n, dtype=np.complex128)
        )

    def to_matrix(self):
        """
        Convert the Multidiags instance to an actual matrix.

        Returns:
            Matrix:
                The matrix representation of the Multidiags instance.
        """
        A = Matrix(CC, self.n, self.n, sparse=True)

        for j, diag in self.diags.items():
            for i, val in enumerate(diag):
                A[i, i + j - self.n] = val
        return A

    def _check_matrix_size(self, other):
        """
        Check if self and other are matrices of the same size.

        Args:
            other (Multidiags):
                Another Multidiags instance.

        Raises:
            TypeError:
                If other is not a Multidiags instance.
            ValueError:
                If matrices are of different sizes.
        """
        if not isinstance(other, Multidiags):
            raise TypeError(
                "The other object must be an instance of Multidiags."
            )
        if self.n != other.n:
            raise ValueError("The matrices must have the same size.")

    def _check_vector_size(self, other):
        """
        Check if other is a vector (numpy array) of correct length.

        Args:
            other (np.ndarray):
                A vector to check.

        Raises:
            TypeError:
                If other is not a numpy array.
            ValueError:
                If the vector has incorrect length.
        """
        if not isinstance(other, np.ndarray):
            raise TypeError("The other object must be a numpy array.")
        if self.n != len(other):
            raise ValueError("The vector is not of correct length.")

    def __add__(self, other):
        """
        Add two Multidiags matrices.

        Args:
            other (Multidiags):
                The matrix to add.

        Returns:
            Multidiags:
                The result of the addition.

        Raises:
            TypeError:
                If other is not a Multidiags instance.
            ValueError:
                If matrices are of different sizes.
        """
        self._check_matrix_size(other)

        diags = {}

        indices0 = self.get_diag_indices()
        indices1 = other.get_diag_indices()
        indices = list(
            set(indices0 + indices1)
        )  # Union of indices0 and indices1

        for i in indices:
            if i not in indices1:
                diags[i] = self.diags[i]
            elif i not in indices0:
                diags[i] = other.diags[i]
            else:
                diags[i] = self.diags[i] + other.diags[i]

        return Multidiags(self.n, diags)

    def __mul__(self, other):
        """
        Multiply the matrix with another matrix or a vector.

        Args:
            other (Multidiags or np.ndarray):
                The matrix or vector to multiply with.

        Returns:
            Multidiags or np.ndarray:
                The result of the multiplication.
        """
        if isinstance(other, np.ndarray):
            self._check_vector_size(other)
            z = sum(
                self.diags[i] * np.roll(other, -i)
                for i in self.get_diag_indices()
            )
            return z

        self._check_matrix_size(other)
        A = Multidiags(self.n)  # Zero matrix
        for i in self.get_diag_indices():
            for j in other.get_diag_indices():
                diags0 = self.diags[i]
                diags1 = np.roll(other.diags[j], -i)
                B = Multidiags(self.n, {i + j: diags0 * diags1})
                A = A + B
        return A

    def BSGS_mult(self, w):
        """
        Perform fast matrix vector multiplication using the Baby Step
        Giant Step algorithm.

        The non zero diagonal indices of the matrix must form an arithmetic
        progression containing zero.

        Args:
            w (np.ndarray):
                The vector to multiply.

        Returns:
            np.ndarray:
                The result of the multiplication.

        Raises:
            ValueError:
                If the vector has incorrect length.
            ValueError:
                If the matrix is not regular.
        """
        if self.n != len(w):
            raise ValueError("The vector does not have the correct length.")

        U = self.get_symm_diag_indices()
        if 0 not in U:
            raise ValueError("The matrix is not regular.")
        t = len(U)
        if t == 1:
            return self * w
        d = U[1] - U[0]
        u_min = min(U)
        u_max = max(U)

        while t & (t - 1) != 0:
            # Ensure that t is a power of two
            t += 1
        half_log_2_t = log(t, 2) / 2
        k0 = 2 ** floor(half_log_2_t)
        k1 = 2 ** ceil(half_log_2_t)

        rotated_w = {
            d * j + u_min: np.roll(w, -d * j - u_min) for j in range(k1)
        }

        z = np.zeros(self.n, dtype=np.complex128)
        for i in range(k0):
            z0 = np.zeros(self.n, dtype=np.complex128)
            a = d * i * k1
            for j in range(k1):
                b = d * j + u_min
                c = a + b
                if c > u_max:
                    # Zero diagonal
                    continue
                z1 = self.get_diag(c)
                z1 = np.roll(z1, a)
                z1 = z1 * rotated_w[b]
                z0 += z1
            z0 = np.roll(z0, -a)
            z += z0
        return z

    def __repr__(self):
        """
        Return a string representation of the matrix.

        Returns:
            str:
                String representation of the matrix.
        """
        return str(self.to_matrix())


# Roots of unity


def get_roots_of_unity(n):
    """
    Get all n-th roots of unity for even n.

    Args:
        n (int):
            The degree of the roots of unity (must be even).

    Returns:
        np.ndarray:
            Array of roots of unity.

    Raises:
        ValueError:
            If n is not even.
    """
    if n % 2 != 0:
        raise ValueError("We only compute roots of unity for even n.")

    folder = os.path.join("Cache", "Roots of Unity")
    file = os.path.join(folder, f"roots_{n}.npy")

    if os.path.exists(file):
        return np.load(file, allow_pickle=True)

    a = 2 * pi * 1j / n
    exponents = np.arange(n // 2)
    roots = np.exp(exponents * a).astype(np.complex128)
    roots = np.concatenate([roots, -roots])

    pickle.dump(roots, open(file, "wb"))
    return roots


# Matrices involved in (homomorphic) encoding & decoding


def _check_power_two(n):
    """
    Check if n is a power of two.

    Args:
        n (int):
            The number to check.

    Raises:
        ValueError:
            If n is not a power of two or less than two.
    """
    if n & (n - 1) != 0:
        raise ValueError("n must be a power of two.")


def bit_rev(i, num_bits):
    """
    Compute the bit reversed of i.

    Args:
        i (int):
            The integer to reverse.
        num_bits (int):
            The number of bits considered.

    Returns:
        int:
            The bit reversed integer.
    """
    rev = 0
    for _ in range(num_bits):
        rev = (rev << 1) | (i & 1)
        i >>= 1
    return rev


def get_E(n, l, inverse=False):
    """
    Get the matrix E_{n, l} or iE_{n, l} for 2**l < n:
    - The product E_{n, 0} * E_{n, 1} * ... is a variant of the DFT matrix.
      After left multiplying by it, the result is given in bit reversed order.
    - The product ... * iE_{n, 1} * iE_{n, 0} is a variant of the inverse
      DFT matrix. One should only left multiply vectors given in bit reversed
      order by it.

    Args:
        n (int):
            The size parameter (must be a power of two).
        l (int):
            The index of the matrix to compute.
        inverse (bool, optional):
            Whether to compute the inverse matrix. Defaults to False.

    Returns:
        Multidiags:
            The matrix E_{n, l} or iE_{n, l}.

    Raises:
        ValueError:
            If l is too large.
    """
    _check_power_two(n)

    folder = os.path.join("Cache", "Matrices E")
    if inverse == False:
        file = os.path.join(folder, f"E_{n}_{l}.npy")
    else:
        file = os.path.join(folder, f"iE_{n}_{l}.npy")

    if os.path.exists(file):
        return np.load(file, allow_pickle=True)

    half_k = 2**l

    if half_k >= n:
        raise ValueError("l is too large.")

    k = 2 * half_k

    roots = get_roots_of_unity(4 * n)
    half_roots = roots / 2

    diag0 = np.zeros(n, dtype=np.complex128)
    diag1 = np.zeros(n, dtype=np.complex128)
    diag2 = np.zeros(n, dtype=np.complex128)

    log_2_n = log(n, 2)

    indices = [
        (2**l * 5 ** bit_rev(i, int(log_2_n) - 1 - l)) % (4 * n)
        for i in range(n // k)
    ]

    if inverse == False:
        for i in range(n):
            i_div, i_mod = divmod(i, k)
            if i_mod < half_k:
                diag0[i] = roots[0]
                diag1[i] = roots[indices[i_div]]
            else:
                diag0[i] = -roots[indices[i_div]]
                diag2[i] = roots[0]
    else:
        for i in range(n):
            i_div, i_mod = divmod(i, k)
            if i_mod < half_k:
                diag0[i] = half_roots[0]
                diag1[i] = half_roots[0]
            else:
                diag0[i] = -half_roots[-indices[i_div]]
                diag2[i] = half_roots[-indices[i_div]]

    if half_k == n // 2:
        E = Multidiags(n, {0: diag0, half_k: diag1 + diag2})
    else:
        E = Multidiags(n, {0: diag0, half_k: diag1, -half_k: diag2})

    pickle.dump(E, open(file, "wb"))
    return E


def get_F(n, l, inverse=False):
    """
    Get the matrix F_{n, l} or iF_{n, l} for 2**l < n:
    - The product F_{n, 0} * F_{n, 1} * ... is a variant of the DFT matrix.
      One should only left multiply vectors given in bit reversed order by it.
    - The product ... * iF_{n, 1} * iF_{n, 0} is a variant of the inverse
      DFT matrix. After left multiplying by it, the result is given in bit
      reversed order.

    Args:
        n (int):
            The size parameter (must be a power of two).
        l (int):
            The index of the matrix to compute.
        inverse (bool, optional):
            Whether to compute the inverse matrix. Defaults to False.

    Returns:
        Multidiags:
            The matrix F_{n, l} or iF_{n, l}.

    Raises:
        ValueError:
            If l is too large.
    """
    _check_power_two(n)

    folder = os.path.join("Cache", "Matrices F")
    if inverse == False:
        file = os.path.join(folder, f"F_{n}_{l}.npy")
    else:
        file = os.path.join(folder, f"iF_{n}_{l}.npy")

    if os.path.exists(file):
        return np.load(file, allow_pickle=True)

    half_k = 2**l

    if half_k >= n:
        raise ValueError("l is too large.")

    k = 2 * half_k
    n_over_k = n // k
    two_n_over_k = 2 * n_over_k

    roots = get_roots_of_unity(4 * n)
    half_roots = roots / 2

    diag0 = np.zeros(n, dtype=np.complex128)
    diag1 = np.zeros(n, dtype=np.complex128)
    diag2 = np.zeros(n, dtype=np.complex128)

    indices = [(k // 2 * 5**i) % (4 * n) for i in range(n_over_k)]

    if inverse == False:
        for i in range(n):
            i_mod = i % two_n_over_k
            if i_mod < n_over_k:
                diag0[i] = roots[0]
                diag1[i] = roots[indices[i_mod]]
            else:
                diag0[i] = -roots[indices[i_mod - n_over_k]]
                diag2[i] = roots[0]
    else:
        for i in range(n):
            i_mod = i % two_n_over_k
            if i_mod < n_over_k:
                diag0[i] = half_roots[0]
                diag1[i] = half_roots[0]
            else:
                diag0[i] = -half_roots[-indices[i_mod - n_over_k]]
                diag2[i] = half_roots[-indices[i_mod - n_over_k]]

    if l == 0:
        F = Multidiags(n, {0: diag0, n_over_k: diag1 + diag2})
    else:
        F = Multidiags(n, {0: diag0, n_over_k: diag1, -n_over_k: diag2})

    pickle.dump(F, open(file, "wb"))
    return F


def group_matrices(matrices, s, right_to_left=False):
    """
    Group a list of matrices A_0, A_1, A_2, ..., into blocks of size s.

    If right_to_left is False, returns
    [..., A_s * ... * A_{2s - 1}, A_0 * ... * A_{s - 1}].
    Otherwise, returns
    [A_{s - 1} * ... * A_0, A_{2s - 1} * ... * A_s, ...].

    Args:
        matrices (list):
            List of Multidiags matrices to group.
        s (int):
            Block size.
        right_to_left (bool, optional):
            Direction of multiplication. Defaults to False.

    Returns:
        list:
            List of grouped matrices.
    """
    k = len(matrices)
    if k == 0:
        return []

    n = matrices[0].n
    num_blocks = ceil(k / s)
    blocks = [0] * num_blocks

    for L in range(num_blocks):
        l_min = L * s
        l_max = min(l_min + s, k)

        A = Multidiags.get_identity(n)

        if right_to_left == False:
            for i in range(l_min, l_max):
                A = A * matrices[i]
            blocks[num_blocks - 1 - L] = A
        else:
            for i in range(l_min, l_max):
                A = matrices[i] * A
            blocks[L] = A

    return blocks


def get_grouped_E(n, s, inverse=False):
    """
    Group the matrices E_{n, l} or iE_{n, l}:
    - If inverse is False, group matrices E_{n, l} as
      [..., E_{n, s} * ... * E_{n, 2s - 1}, E_{n, 0} * ... * E_{n, s - 1}].
    - Otherwise, group matrices iE_{n, l} as
      [iE_{n, s - 1} * ... * iE_{n, 0}, iE_{n, 2s - 1} * ... * iE_{n, s}, ...].

    Args:
        n (int):
            The size parameter (must be a power of two).
        s (int):
            Block size.
        inverse (bool, optional):
            Whether to compute the inverse matrices. Defaults to False.

    Returns:
        list:
            List of grouped matrices.
    """
    _check_power_two(n)

    folder = os.path.join("Cache", "Grouped Matrices E")
    if inverse == False:
        file = os.path.join(folder, f"grouped_E_{n}_{s}.npy")
    else:
        file = os.path.join(folder, f"grouped_iE_{n}_{s}.npy")

    if os.path.exists(file):
        return np.load(file, allow_pickle=True)

    matrices = [get_E(n, l, inverse) for l in range(log(n, 2))]
    grouped_matrices = group_matrices(matrices, s, inverse)

    pickle.dump(grouped_matrices, open(file, "wb"))
    return grouped_matrices


def get_grouped_F(n, s, inverse=False):
    """
    Group the matrices F_{n, l} or iF_{n, l}:
    - If inverse is False, group matrices F_{n, l} as
      [..., F_{n, s} * ... * F_{n, 2s - 1}, F_{n, 0} * ... * F_{n, s - 1}].
    - Otherwise, group matrices iF_{n, l} as
      [iF_{n, s - 1} * ... * iF_{n, 0}, iF_{n, 2s - 1} * ... * iF_{n, s}, ...].

    Args:
        n (int):
            The size parameter (must be a power of two).
        s (int):
            Block size.
        inverse (bool, optional):
            Whether to compute the inverse matrices. Defaults to False.

    Returns:
        list:
            List of grouped matrices.
    """
    _check_power_two(n)

    folder = os.path.join("Cache", "Grouped Matrices F")
    if inverse == False:
        file = os.path.join(folder, f"grouped_F_{n}_{s}.npy")
    else:
        file = os.path.join(folder, f"grouped_iF_{n}_{s}.npy")

    if os.path.exists(file):
        return np.load(file, allow_pickle=True)

    matrices = [get_F(n, l, inverse) for l in range(int(log(n, 2)))]
    grouped_matrices = group_matrices(matrices, s, inverse)

    pickle.dump(grouped_matrices, open(file, "wb"))
    return grouped_matrices
