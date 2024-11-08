
# CKKS in SageMath
A CKKS implementation in SageMath which supports basic bootstrapping.

This repository provides a SageMath implementation of the CKKS homomorphic encryption scheme [2], designed for efficient approximate computations on encrypted data. CKKS is particularly suited for applications involving complex numbers and approximate arithmetic.

## Table of Contents
- [Overview](#overview)
- [Installation](#installation)
- [Files and Structure](#files-and-structure)
- [Usage](#usage)
- [Class Overview](#class-overview)
- [Contributing](#contributing)
- [License](#license)
- [References](#references)

## Overview

The CKKS scheme enables approximate arithmetic operations on encrypted data, which is useful for privacy preserving computations, where exact results are unnecessary or impossible to achieve. This implementation is based on SageMath. The scheme was originally introduced by Cheon et al. in 2016 [1], with subsequent work on bootstrapping by Cheon et al. in 2018 [3] and running time improvements by Cheon et al. in 2018 [4].

## Installation

1. **Prerequisites**: Ensure you have SageMath 10.0 or higher installed. You can download it from the [SageMath website](https://www.sagemath.org/download.html).
   
2. **Clone the Repository**:
   ```bash
   git clone https://github.com/se-tim/CKKS-in-SageMath.git
   cd CKKS-in-SageMath
   ```

3. **Dependencies**: There are no additional dependencies beyond SageMath.

## Files and Structure

- **ckks.py**: Contains the CKKS class, implementing key aspects of the CKKS scheme, including homomorphic arithmetic and bootstrapping [3].
  
- **poly.py**: Manages polynomial operations required for CKKS, like polynomial multiplication and modular reduction.

- **fast_dft.py**: Implements a fast Discrete Fourier Transform (DFT), optimizing polynomial multiplications within the CKKS scheme [4].

- **test.ipynb**: A Jupyter notebook demonstrating CKKS functionality with examples for setting parameters, key generation, encryption, decryption, homomorphic operations, and bootstrapping.

- **Cache**: This directory contains cached files generated during runtime. Do not delete.

- **Estimator**: Security estimation for RLWE and LWE instances [1].

## Usage

### Running Tests and Examples

The `test.ipynb` notebook demonstrates the following steps in CKKS encryption:

1. **Configuration**: Set CKKS parameters.
2. **Key Generation**: Generate keys for encryption and decryption.
3. **Encoding and Encryption**: Encode complex vectors into polynomials and encrypt them.
4. **Decryption and Decoding**: Verify encryption by decrypting and decoding ciphertexts.
5. **Homomorphic Operations**: Perform addition and multiplication on ciphertexts and verify results.
6. **Bootstrapping**: Refresh ciphertexts to raise levels, and verify with decryption.

Run the notebook cells sequentially to see these operations in action.

## Class Overview

### CKKS Class Functionalities

- **Configuration (`config`)**: Set key CKKS parameters such as ring degree, slot count, maximum level, smallest modulus, and scaling factor.
  
- **Key Generation (`key_gen`)**: Generates the secret, public, and evaluation keys, with a security level estimation.

- **Encoding/Decoding**: 
  - **`encode`**: Encodes complex vectors as integer polynomials.
  - **`decode`**: Decodes integer polynomials back into complex vectors.

- **Encryption/Decryption**:
  - **`enc_poly_with_sk`**: Encrypts a plaintext polynomial using the secret key.
  - **`dec_to_poly`**: Decrypts ciphertext back to a plaintext polynomial.

- **Homomorphic Operations**:
  - **Arithmetic** (`__add__`, `__mul__`): Supports addition, multiplication, and other basic arithmetic between ciphertexts or between ciphertexts and plaintexts.
  - **Modular Reduction and Rescaling** (`__mod__`, `rescale`): Allows modulus switching and rescaling to control noise.
  - **Rotation, Conjugation, and Galois Automorphisms** (`rotate`, `conjugate`, `galois`): Enables Galois automorphisms and vector rotations, which are essential for operations on encrypted vectors.

- **Bootstrapping**:
  - **`config_bootstrap`**: Precomputes values required for bootstrapping.
  - **`bootstrap`**: Refreshes ciphertexts to raise their levels.

## License

This project is open source and available under the MIT License. See `LICENSE` for more details.

## References

1. Albrecht, M.R., Player, R., & Scott, S. (2015). *On the Concrete Hardness of Learning with Errors*. Journal of Mathematical Cryptology, Volume 9, Issue 3, Pages 169–203. [doi:10.1515/jmc-2015-0016](https://doi.org/10.1515/jmc-2015-0016)

2. Cheon, J.H., Kim, A., Kim, M., & Song, Y. (2016). *Homomorphic Encryption for Arithmetic of Approximate Numbers*. Cryptology ePrint Archive, Paper 2016/421. [https://eprint.iacr.org/2016/421](https://eprint.iacr.org/2016/421)

3. Cheon, J.H., Han, K., Kim, A., Kim, M., & Song, Y. (2018). *Bootstrapping for Approximate Homomorphic Encryption*. Cryptology ePrint Archive, Paper 2018/153. [https://eprint.iacr.org/2018/153](https://eprint.iacr.org/2018/153)

4. Cheon, J.H., Han, K., & Hhan, M. (2018). *Faster Homomorphic Discrete Fourier Transforms and Improved FHE Bootstrapping*. Cryptology ePrint Archive, Paper 2018/1073. [https://eprint.iacr.org/2018/1073](https://eprint.iacr.org/2018/1073)
