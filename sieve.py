from math import ceil


check_quadres = lambda n, p: p == 2 or pow(n, (p - 1) // 2, p) == 1


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


def quadsieve(n, B, k):
    '''Find B-smooth numbers and their exponent vectors out of
    the sequence x^2 - n for x from sqrt(n) to sqrt(n)+k.
    Return a tuple (x, x^2 - n, v) where v is the exponent vector.
    '''
    primes = sieve(B)
    filtered_primes = [p for p in primes if check_quadres(n % p, p)]
    minx = ceil(n**0.5)
    bsize = len(filtered_primes)
    seqx = [x + minx for x in range(k)]
    seqq = [x ** 2 - n for x in seqx]
    max_x = minx + k - 1
    seqf = [bsize * [0] for _ in range(k)]
    for ip, p in enumerate(filtered_primes):
        sol = solve_quadres(n % p, p)
        sc1 = sol + p * ceil((minx - sol) / p)
        for s in range(sc1, max_x + 1, p):
            i = s - minx
            while not (seqq[i] % p):
                seqq[i] //= p
                seqf[i][ip] = (seqf[i][ip] + 1) % 2  
        sol2 = p - sol
        sc2 = sol2 + p * ceil((minx - sol2) / p)
        for s in range(sc2, max_x + 1, p):
            i = s - minx
            while not (seqq[i] % p):
                seqq[i] //= p
                seqf[i][ip] = (seqf[i][ip] + 1) % 2
    return [(x, x**2 - n, v,) for x, t, v in zip(seqx, seqq, seqf) if t == 1]
