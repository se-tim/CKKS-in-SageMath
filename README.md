
# CKKS in SageMath
A CKKS implementation in SageMath which supports basic bootstrapping.

This repository provides a SageMath implementation of the CKKS homomorphic encryption scheme [2],
designed for efficient approximate computations on encrypted data.
CKKS is particularly suited for applications involving complex numbers and approximate arithmetic.

## Table of Contents
- [Overview](#overview)
- [Installation](#installation)
- [Repository Structure](#files-and-structure)
- [Usage](#usage)
- [License](#license)
- [References](#references)

## Overview

The CKKS scheme enables approximate arithmetic operations on encrypted data,
which is useful for privacy preserving computations
for which exact results are unnecessary or impossible to achieve.
The scheme was originally introduced by Cheon et al. in 2016 [2],
with subsequent work on bootstrapping by Cheon et al. in 2018 [3]
and running time improvements for bootstrapping by Cheon et al. in 2018 [4].
This SageMath implementation is based on these papers.

## Installation

1. **Prerequisites**: Ensure you have SageMath 10.0 or higher installed.
You can download it from the [SageMath website](https://www.sagemath.org/download.html).
   
2. **Clone the repository**:
   ```bash
   git clone https://github.com/se-tim/CKKS-in-SageMath.git
   cd CKKS-in-SageMath
   ```

3. **Initialize the submodules**: If you intend to use this CKKS implementation as a package,
you must ensure all submodules are correctly initialized,
as explained here:
[Git submodule inside of a submodule (nested submodules)](https://stackoverflow.com/q/1535524).
Therefore, run the following command inside the cloned repository:
   ```bash
   git submodule update --init --recursive
   ```

## Repository Structure

- **`ckks_package/`**: The main package containing all implementations related to the CKKS scheme.

  - **`bit_rev.py`**: Implements functions for bit reversal operations.

  - **`ckks.py`**: Contains the CKKS class,
  implementing key aspects of the CKKS scheme,
  including homomorphic arithmetic [2] and bootstrapping [3].

  - **`fast_dft.py`**: Implements a version of the Discrete Fourier Transform [4],
  essential for encoding and decoding complex vectors in CKKS.
  It is also required for the CoeffToSlot and SlotToCoeff transformations,
  both of which are necessary for bootstrapping in CKKS.

  - **`poly.py`**: Manages polynomial operations required for CKKS.

  - **`lattice_estimator/`**: This is a submodule.
  It is a tool for security estimation in the context of RLWE and LWE [1].
  Make sure that this is initialized correctly,
  as described in the section [Installation](#installation).

- **`test.ipynb`**: A Jupyter notebook demonstrating the CKKS functionality with examples.  

## Usage

The `test.ipynb` notebook demonstrates the key components of the CKKS scheme, including:

1. **Configuration**: Set CKKS parameters.
2. **Key Generation**: Generate keys for encryption and decryption.
3. **Encoding and Encryption**: Encode complex vectors into polynomials and encrypt them.
4. **Decryption and Decoding**: Verify encryption by decrypting and decoding ciphertexts.
5. **Homomorphic Operations**: Perform addition and multiplication on ciphertexts,
and verify the results.
6. **Bootstrapping**: Refresh a ciphertext to raise its level,
and verify the result.

Run the notebook cells sequentially to see these operations in action.

## License

This project is open source and available under the MIT License.
See `LICENSE` for more details.

## References

1. Martin R. Albrecht, Rachel Player, & Sam Scott.
*On the Concrete Hardness of Learning with Errors*.
Journal of Mathematical Cryptology, Volume 9, Issue 3, Pages 169–203, 2015.
[doi:10.1515/jmc-2015-0016](https://doi.org/10.1515/jmc-2015-0016)

2. Jung Hee Cheon, Andrey Kim, Miran Kim, & Yongsoo Song.
*Homomorphic Encryption for Arithmetic of Approximate Numbers*.
Cryptology ePrint Archive, Paper 2016/421, 2016.
[https://eprint.iacr.org/2016/421](https://eprint.iacr.org/2016/421)

3. Jung Hee Cheon, Kyoohyung Han, Andrey Kim, Miran Kim, & Yongsoo Song.
*Bootstrapping for Approximate Homomorphic Encryption*.
Cryptology ePrint Archive, Paper 2018/153, 2018.
[https://eprint.iacr.org/2018/153](https://eprint.iacr.org/2018/153)

4. Jung Hee Cheon, Kyoohyung Han, & Minki Hhan.
*Faster Homomorphic Discrete Fourier Transforms and Improved FHE Bootstrapping*.
Cryptology ePrint Archive, Paper 2018/1073, 2018.
[https://eprint.iacr.org/2018/1073](https://eprint.iacr.org/2018/1073)
