from math import ceil
from quadres import solve_quadres, check_quadres


def sieve(n):
    '''Return a list of primes up to and including n.
    Implements the sieve of Eratosthenes.
    '''
    A = [True] * (n - 1)
    for i in range(2, ceil(n**0.5)):
        if A[i - 2]:
            A[2 * i - 2::i] = [False] * len(A[2 * i - 2::i])
    return [i + 2 for i in range(len(A)) if A[i]]


def quadsieve(n, B, k, add_terms=False):
    '''Find B-smooth numbers and their exponent vectors out of
    the sequence x^2 - n for x from sqrt(n) to sqrt(n)+k.
    Return a tuple (x, x^2 - n, v) where v is the exponent vector.
    '''
    primes = sieve(B)
    filtered_primes = [p for p in primes if check_quadres(n % p, p)]
    x_centre = ceil(n**0.5)
    # pi(B)
    bsize = len(filtered_primes)
    seqx = [x for x in range(max(x_centre - k, 1), x_centre + k)]
    seqq = [x ** 2 - n for x in seqx]
    max_x = seqx[-1]
    min_x = seqx[0]
    # Last entry is for the sign of the number
    seqf = [bsize * [0] + [0 if x >= 0 else 1] for x in seqq]
    gens = []
    for ip, p in enumerate(filtered_primes):
        # Solve x^2 = n mod p
        # We only get the first solution here
        sol = solve_quadres(n % p, p)
        # Smallest value of x >= min_x such that x = a mod p
        sc1 = sol + p * ceil((min_x - sol) / p)
        for s in range(sc1, max_x + 1, p):
            i = s - min_x
            # Repeatedly divide x^2 - n by p until it is no longer divisible
            while not (seqq[i] % p):
                seqq[i] //= p
                seqf[i][ip] = (seqf[i][ip] + 1) % 2

        # Repeat the steps above for the other solution
        sol2 = p - sol
        sc2 = sol2 + p * ceil((min_x - sol2) / p)
        for s in range(sc2, max_x + 1, p):
            i = s - min_x
            while not (seqq[i] % p):
                seqq[i] //= p
                seqf[i][ip] = (seqf[i][ip] + 1) % 2
        # Add these solutions to an array for later use
        # if we need to find more B-smooth terms.
        if add_terms:
            gens.append((p, ip, sol,))
            gens.append((p, ip, sol2,))
    vectors = [(x, x**2 - n, v) for x, t, v in
               zip(seqx, seqq, seqf) if abs(t) == 1]
    # Find enough B-smooth terms such that linear dependence is guaranteed.
    # Only runs if the add_terms is set to True.
    k = max_x + 1
    kn = min_x - 1
    if add_terms:
        while len(vectors) < bsize + 2:
            q = k ** 2 - n
            mq = q
            v = bsize * [0] + [0]
            for p, ip, sol in gens:
                while not (mq % p):
                    mq //= p
                    v[ip] = (v[ip] + 1) % 2
            if abs(mq) == 1:
                vectors.append((k, q, v,))
            k += 1

            # We don't care about negative x since (-x)^2 = x^2
            if kn >= 1:
                q = kn ** 2 - n
                mq = q
                v = bsize * [0] + [1]
                for p, ip, sol in gens:
                    while not (mq % p):
                        mq //= p
                        v[ip] = (v[ip] + 1) % 2
                if abs(mq) == 1:
                    vectors.append((kn, q, v,))
                kn -= 1

    return vectors
