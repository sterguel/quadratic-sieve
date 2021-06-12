from math import ceil


check_quadres = lambda n, p: p == 2 or pow(n, (p - 1) // 2, p) == 1
gcd = lambda a, b: gcd(b, a % b) if b != 0 else a


def solve_quadres(n, p):
    '''Return an integer r such that r^2 = n mod p (if it exists).
    Implements the Tonelli-Shanks algorithm.
    '''
    if not check_quadres(n, p):
        raise ValueError(f'{n} has no square root in Z/{p}Z')
    if p == 2:
        return n
    # Find Q,S such that p - 1 = Q*2^S
    Q = p - 1
    S = 0
    while not (Q % 2):
        Q //= 2
        S += 1
    z = 2
    while 1:
        if not check_quadres(z, p):
            break
        else:
            z += 1
    M = S
    c = pow(z, Q, p)
    t = pow(n, Q, p)
    R = pow(n, (Q + 1) // 2, p)
    while True:
        if not t:
            return 0
        elif t == 1:
            return R
        mi = 0
        for i in range(1, M):
            if pow(t, 2**i, p) == 1:
                mi = i
                break
        eb = pow(2, M - mi - 1)
        b = pow(c, eb, p)
        M = i
        c = pow(b, 2, p)
        t = (t * c) % p
        R = (R * b) % p


def sieve(n):
    '''Return a list of primes up to and including n'''
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
    bsize = len(filtered_primes)
    seqx = [x for x in range(max(x_centre - k, 1), x_centre + k)]
    seqq = [x ** 2 - n for x in seqx]
    max_x = seqx[-1]
    min_x = seqx[0]
    # First entry is for the sign of the number
    seqf = [[0 if x >= 0 else 1] + bsize * [0] for x in seqq]
    gens = []
    for ip, p in enumerate(filtered_primes):
        sol = solve_quadres(n % p, p)
        sc1 = sol + p * ceil((min_x - sol) / p)
        for s in range(sc1, max_x + 1, p):
            i = s - min_x
            while not (seqq[i] % p):
                seqq[i] //= p
                seqf[i][ip + 1] = (seqf[i][ip + 1] + 1) % 2
        sol2 = p - sol
        sc2 = sol2 + p * ceil((min_x - sol2) / p)
        for s in range(sc2, max_x + 1, p):
            i = s - min_x
            while not (seqq[i] % p):
                seqq[i] //= p
                seqf[i][ip + 1] = (seqf[i][ip + 1] + 1) % 2
        if add_terms:
            gens.append((p, ip, sol,))
            gens.append((p, ip, sol2,))
    print(seqq)
    vectors = [(x, x**2 - n, v) for x, t, v in
               zip(seqx, seqq, seqf) if abs(t) == 1]
    k = max_x + 1
    if add_terms:
        while len(vectors) < bsize + 2:
            q = k ** 2 - n
            mq = q
            v = [0] + bsize * [0]
            for p, ip, sol in gens:
                while not (mq - sol) % p:
                    mq //= p
                    v[ip + 1] = (v[ip + 1] + 1) % 2
            if abs(mq) == 1:
                vectors.append((k, q, v,))
            k += 1
    return vectors
