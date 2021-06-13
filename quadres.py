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
