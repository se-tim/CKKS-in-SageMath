
# CKKS in SageMath

This repository provides a SageMath implementation
of the CKKS homomorphic encryption scheme [2],
designed for efficient approximate computations on encrypted data.
CKKS is particularly suited for applications
involving complex numbers and approximate arithmetic.

## Table of Contents

- [Installation](#installation)
- [Repository Structure](#repository-structure)
- [Usage](#usage)
- [License](#license)
- [References](#references)

## Installation

1. **Prerequisites**:
Ensure you have SageMath 10.0 or higher installed.
You can download it from the
[SageMath website](https://www.sagemath.org/download.html).
   
2. **Clone the repository**:
   ```bash
   git clone https://github.com/se-tim/CKKS-in-SageMath.git
   cd CKKS-in-SageMath
   ```

3. **Initialize the submodules**:
If you intend to use this CKKS implementation as a package,
you must ensure all submodules are correctly initialized,
as explained here:
[Git submodule inside of a submodule (nested submodules)](https://stackoverflow.com/q/1535524).
Therefore, run the following command inside the cloned repository:
   ```bash
   git submodule update --init --recursive
   ```

## Repository Structure

- **`ckks_package/`**:
The main package containing all implementations related to the CKKS scheme:

  - **`bit_rev.py`**:
  Implements functions for bit-reversal operations.

  - **`ckks.py`**:
  Contains the CKKS class,
  implementing key aspects of the CKKS scheme,
  including homomorphic arithmetic [2] and bootstrapping [3].

  - **`fast_dft.py`**:
  Implements a version of the Discrete Fourier Transform [4],
  essential for encoding and decoding complex vectors in CKKS.
  It is also required for the CoeffToSlot and SlotToCoeff transformations,
  both of which are necessary for bootstrapping in CKKS.

  - **`poly.py`**:
  Manages polynomial operations required for CKKS.

  - **`lattice_estimator/`**:
  This is a submodule.
  It is a tool for security estimation in the context of RLWE and LWE [1].
  Make sure that this is initialized correctly,
  as described in the section [Installation](#installation).

- **`test.py`**
 A script demonstrating the main CKKS functionality.

## Usage

The `test.py` script demonstrates the key components of the CKKS scheme
through a linear sequence of operations:

1. **Parameter configuration**:
Choose default parameters or manually input your own
for ring degree, number of slots, maximal level,
base modulus, and scaling factors.

2. **Key generation**:
Generate secret, public and evaluation keys.

3. **Security estimation**:
Estimate the security level of the current parameter set.
By default, the primal hybrid security level is not estimated (to save time).

4. **Encryption**:
Encrypt random complex vectors.

5. **Homomorphic operations**:
Perform homomorphic addition and multiplication on the ciphertexts.

6. **Bootstrapping**:
Refresh a ciphertext at lowest level through bootstrapping.

7. **Precision measurement**:
Evaluate the precision loss introduced during bootstrapping.

To run the script,
make sure you are in the root directory of the repository,
then execute:
```bash
sage test.py
```
It will output step-by-step progress and results directly to the terminal.

## License

This project is open source and available under the MIT License.
See `LICENSE` for more details.

## References

1. Martin Albrecht, Rachel Player, and Sam Scott.
On the Concrete Hardness of Learning with Errors.
*Journal of Mathematical Cryptology*, 9, 2015.

2. Jung Hee Cheon, Dongwoo Kim, Duhyeong Kim, and Yongsoo Song.
Homomorphic Encryption for Approximate Numbers.
In *Advances in Cryptology – ASIACRYPT 2017*.
Springer, 2017.

3. Jung Hee Cheon, Kyoohyung Han, Andrey Kim, Miran Kim, and Yongsoo Song.
Bootstrapping for Approximate Homomorphic Encryption.
In *Advances in Cryptology – EUROCRYPT 2018*.
Springer, 2018.

4. Jung Hee Cheon, Kyoohyung Han, and Minki Hhan.
Improved Homomorphic Discrete Fourier Transforms and FHE Bootstrapping.
*IEEE Access*, 7, 2019.