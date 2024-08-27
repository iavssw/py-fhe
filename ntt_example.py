

from math import pi, cos, sin
import os
import unittest
from util.ntt import NTTContext, FFTContext
from util.polynomial import Polynomial
from util.random_sample import sample_uniform
from tests.helper import check_complex_vector_approx_eq

def multiply_naive(a, b, ring_degree, coeff_modulus=None):

    poly_prod = Polynomial(ring_degree, [0] * ring_degree)

    for d in range(2 * ring_degree - 1):
        # Since x^d = -1, the degree is taken mod d, and the sign
        # changes when the exponent is > d.
        index = d % ring_degree
        sign = int(d < ring_degree) * 2 - 1

        # Perform a convolution to compute the coefficient for x^d.
        coeff = 0
        for i in range(ring_degree):
            if 0 <= d - i < ring_degree:
                coeff += a[i] * b[d - i]
        poly_prod.coeffs[index] += sign * coeff
        
        if coeff_modulus:
            poly_prod.coeffs[index] %= coeff_modulus

    return poly_prod


# setup parameters

ring_degree = 4 # assume for FHE this would be 1024
coeff_modulus_prime = 7681 # assume for FHE this would be 2^512 prime

coef1 = [1, 2, 3, 4]
coef2 = [5, 6, 7, 8]


# Naive Ring Polynomial multiplication

naive_mul = multiply_naive(coef1, coef2, ring_degree,  coeff_modulus=coeff_modulus_prime)
print("Naive multiplication: ", naive_mul)


# NTT polynomial multiplication

ntt = NTTContext(poly_degree=ring_degree, coeff_modulus=coeff_modulus_prime)

a = ntt.ftt_fwd(coef1)
b = ntt.ftt_fwd(coef2)
ab = [a[i] * b[i] for i in range(ring_degree)]
prod = ntt.ftt_inv(ab)
prod = Polynomial(ring_degree, prod)

print("NTT multiplication: ", prod)