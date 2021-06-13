from sieve import quadsieve
from linalg import find_linear_dependence
from quadres import gcd
from math import isqrt


def factor(N, B, M, add_vectors=True):
    vectors = quadsieve(N, B, M, add_terms=add_vectors)
    nullspace = find_linear_dependence([v[2] for v in vectors])
    if not nullspace:
        raise ValueError('Not enough exponent vectors found.\
Try increasing B or M.')
    divisors = set()
    for b in nullspace:
        x = 1
        y = 1
        for v, e in zip(vectors, b):
            if e:
                x = (x * v[0]) % N
                y *= v[1]
        y = isqrt(y) % N
        divisors.add(abs(gcd(x - y, N)))
        divisors.add(abs(gcd(x + y, N)))
    return divisors - {1, N}
