from sage.all import log, randint
from ckks_package.ckks import CKKS
from ckks_package.poly import Poly

# Configuration

if input("Want to enter your own parameters? (y/n) ") == "y":
    N = eval(input("Ring degree: N = "))
    n = eval(input("Number of slots: n = "))
    L_boot = eval(input("Maximal level during bootstrapping: L_boot = "))
    q = eval(input("Base modulus: q = "))
    p = eval(input("Scaling factor outside of bootstrapping: p = "))
    delta = eval(input("Scaling factor during bootstrapping: delta = "))
else:
    N = 2**15
    n = 2**2
    L_boot = 20
    q = 2**28
    p = 2**22
    delta = 2**32
    print("We chose the following parameters:")
    print(f"N = 2^{log(N,2)} (ring degree)")
    print(f"n = 2^{log(n,2)} (number of slots)")
    print(f"L_boot = {L_boot} (maximal level during bootstrapping)")
    print(f"q = 2^{log(q,2)} (base modulus)")
    print(f"p = 2^{log(p,2)} (scaling factor outside of bootstrapping)")
    print(f"delta = 2^{log(delta,2)} (scaling factor during bootstrapping)")
print()

CKKS.config(N, n, L_boot, q, p, delta, print_messages=True)
print()

CKKS.key_gen(print_messages=True)
print()

CKKS.config_bootstrap(CKKS.sk, print_messages=True)
print()

print("Estimated security in bits:")
print(CKKS.get_security(check_primal_hybrid=False))
print()

# Testing main functionalities

print("Two randomly generated complex vectors:")
polynomials = [
    Poly(
        [randint(-p, p) if i % (N // (2 * n)) == 0 else 0 for i in range(N)],
        N,
    )
    for _ in range(2)
]
vectors = [CKKS.decode(polynomials[i]) for i in range(2)]
print(vectors[0])
print(vectors[1])
print()

print("Encrypting them:")
polynomials = [CKKS.encode(z) for z in vectors]
ciphertexts = [CKKS.enc_poly_with_pk(pol, CKKS.pk) for pol in polynomials]
print(ciphertexts[0])
print(ciphertexts[1])
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
