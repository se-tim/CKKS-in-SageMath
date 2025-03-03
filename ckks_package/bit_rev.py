import numpy as np
from sage.all import *


def bit_rev(i, num_bits):
    """
    Bit reverse i.

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


def bit_rev_vector(z, num_bits=None):
    """
    Rearrange the entries of a vector into bit reversed order.

    Args:
        z (np.ndarray):
            The input vector whose entries should be rearranged.
        num_bits (int, optional):
            The number of bits considered for bit reversal. Defaults to
            log(len(z), 2).

    Returns:
        np.ndarray:
            A new vector with entries rearranged in bit reversed order.
    """
    if num_bits is None:
        num_bits = log(len(z), 2)
    return np.array([z[bit_rev(i, num_bits)] for i in range(len(z))])
