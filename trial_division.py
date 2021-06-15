from math import isqrt


def trial_division(N):
    '''Find a factor of N by dividing by consecutive
    prime numbers until one is found.
    '''
    sq = isqrt(N)
    sieve = [True] * (sq - 1)
    for k in range(2, sq + 1):
        i = k - 2
        if sieve[i]:
            sieve[i + k::k] = [False] * len(sieve[i + k::k])
            if not (N % k):
                return k


def trial_division_bitpacked(N):
    '''Find a factor of N by dividing by consecutive
    prime numbers until one is found.
    Uses bit packing to avoid running out of memory.
    '''
    sq = isqrt(N)
    sieve = (1 << sq) - 1
    for k in range(2, sq + 1):
        i = k - 2
        if sieve & (1 << i):
            sieve &= ~sum([1 << j for j in range(i + k, sq, k)])
            if not (N % k):
                return k
