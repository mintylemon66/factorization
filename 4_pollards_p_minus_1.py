import math
import random
import time

# Pollard's p-1, single stage.
# Parameters: B = 100,000; bases 2, 3, 5.
# M is built from all prime powers q^e <= B (so the method finds a prime
# factor p whenever p-1 is B-powersmooth: every prime power in p-1 is <= B).
# gcd is checked every 200 primes; on overshoot (gcd == n) we backtrack from
# the last checkpoint one multiplication at a time, which separates factors
# unless every factor surfaces at the same single multiplication.
# Added July 2026, alongside the Pollard's rho implementation used for the
# original timings (see 3_pollards_rho.py and the README).

B = 100_000

def _sieve(limit):
    bs = [True] * (limit + 1)
    bs[0] = bs[1] = False
    for p in range(2, int(limit ** 0.5) + 1):
        if bs[p]:
            bs[p * p::p] = [False] * len(bs[p * p::p])
    return [i for i, b in enumerate(bs) if b]

_PRIMES = _sieve(B)
_PP = [(q, q ** int(math.log(B) / math.log(q))) for q in _PRIMES]

def is_probable_prime(a: int, k: int = 10) -> bool:
    if a < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if a % p == 0:
            return a == p
    d, r = a - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for _ in range(k):
        x = pow(random.randrange(2, a - 1), d, a)
        if x in (1, a - 1):
            continue
        for _ in range(r - 1):
            x = x * x % a
            if x == a - 1:
                break
        else:
            return False
    return True

def p_minus_1(n: int, bases=(2, 3, 5), check_every: int = 200):
    """Return a nontrivial factor of composite n, or None on failure."""
    if n % 2 == 0:
        return 2
    for b in bases:
        g = math.gcd(b, n)
        if 1 < g < n:
            return g
        if g == n:
            continue
        x = b
        ckpt, ckpt_i = x, 0
        overshoot = False
        for idx in range(0, len(_PP), check_every):
            for q, qe in _PP[idx:idx + check_every]:
                x = pow(x, qe, n)
            d = math.gcd(x - 1, n)
            if d == n:
                # backtrack from checkpoint, one multiplication at a time
                y = ckpt
                stop = False
                for q, qe in _PP[ckpt_i:idx + check_every]:
                    e = 0
                    while q ** (e + 1) <= qe:
                        e += 1
                    for _ in range(e + 1):
                        y = pow(y, q, n)
                        d2 = math.gcd(y - 1, n)
                        if 1 < d2 < n:
                            return d2
                        if d2 == n:
                            stop = True
                            break
                    if stop:
                        break
                overshoot = True
                break
            if d > 1:
                return d
            ckpt, ckpt_i = x, idx + check_every
        if overshoot:
            continue
        d = math.gcd(x - 1, n)
        if 1 < d < n:
            return d
    return None

if __name__ == "__main__":
    # timing harness, same convention as the other scripts:
    # paste a dataset as my_list = [...] and run; prints one CPU-ms per input.
    my_list = [2 ** i - 1 for i in range(2, 108)]
    for N in my_list:
        start_time = time.process_time()
        result = N if is_probable_prime(N) else p_minus_1(N)
        elapsed_ms = (time.process_time() - start_time) * 1000
        print(f"{elapsed_ms}\t{'ok' if result is not None else 'FAIL'}")
