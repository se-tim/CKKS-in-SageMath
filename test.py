from sage.all import log
from ckks_package.ckks import CKKS


# Configuration

if input("Enter your own parameters (y / n)? ") == "y":
    N = 2 ** int(input("Ring degree: N = 2^"))
    n = 2 ** int(input("Number of slots: n = 2^"))
    L = int(input("Maximal level during bootstrapping: L = "))
    q = 2 ** int(input("Base modulus: q = 2^"))
    p = 2 ** int(input("Scaling factor outside of bootstrapping: p = 2^"))
    delta = 2 ** int(input("Scaling factor during bootstrapping: delta = 2^"))
else:
    N = 2**15
    n = 2**2
    L = 20
    q = 2**28
    p = 2**22
    delta = 2**32
    print("We chose the following parameters:")
    print(f"Ring degree: N = 2^{log(N, 2)}")
    print(f"Number of slots: n = 2^{log(n, 2)}")
    print(f"Maximal level during bootstrapping: L = {L}")
    print(f"Base modulus: q = 2^{log(q, 2)}")
    print(f"Scaling factor outside of bootstrapping: p = 2^{log(p, 2)}")
    print(f"Scaling factor during bootstrapping: delta = 2^{log(delta, 2)}")
print()

CKKS.config(N, n, L, q, p, delta, print_messages=True)
print()

CKKS.key_gen(print_messages=True)
print()

CKKS.config_bootstrap(CKKS.sk, print_messages=True)
print()

print("Estimated security in bits:")
print(CKKS.get_security(check_primal_hybrid=False))
print()

# Testing main functionalities

print("Two randomly generated ciphertexts:")
ciphertexts = [CKKS.get_random_ciphertext(CKKS.sk) for _ in range(2)]
print(ciphertexts[0])
print(ciphertexts[1])
print()

print("They decrypt to:")
for i in range(2):
    print(CKKS.decode(ciphertexts[i].dec_to_poly(CKKS.sk)))
print()

print("Homomorphic addition and multiplication:")
ct_add = ciphertexts[0] + ciphertexts[1]
ct_mult = ciphertexts[0] @ ciphertexts[1]
print(ct_add)
print(ct_mult)
print()

print("They decrypt to:")
print(CKKS.decode(ct_add.dec_to_poly(CKKS.sk)))
print(CKKS.decode(ct_mult.dec_to_poly(CKKS.sk)))
print()

print("Projecting down the addition ciphertext to the lowest level:")
ct_low = ct_add % q
print(ct_low)
print()

print("Bootstrapping it:")
ct_boot = ct_low.bootstrap()
print(ct_boot)
print()

print("This decrypts to:")
print(CKKS.decode(ct_boot.dec_to_poly(CKKS.sk)))
print()

print("Precision in bits:")
print(ct_low.get_precision(ct_boot, CKKS.sk))
