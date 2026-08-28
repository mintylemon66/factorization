# Base-failure density on 2047 = 23 x 89. For every valid base b (2 <= b <= 2045,
# gcd(b, 2047) = 1) this measures whether the staged p-1 run fails (first
# nontrivial gcd is 2047 itself) and compares two exposure-key models: the
# corrected one, where orders 1 and 2 share a key because the first gcd check
# follows the first squaring, and the old one that gave order 1 its own class.
# Expected: 1934 valid bases, 1602 measured failures (0.8283); corrected model
# 1602 with zero disagreements; old model 1600, wrong exactly on b = 622 and 1425.

import math

N = 2047
FACTORS = (23, 89)
B = 100_000


def primes_up_to(limit):
    sieve = [True] * (limit + 1)
    sieve[0] = sieve[1] = False
    for p in range(2, int(limit ** 0.5) + 1):
        if sieve[p]:
            sieve[p * p :: p] = [False] * len(sieve[p * p :: p])
    return [i for i, ok in enumerate(sieve) if ok]


PRIMES = primes_up_to(B)

# The multiplication schedule: for each prime p (ascending), multiply the exponent by p
# one factor at a time, up to the largest p^k <= B.
SCHEDULE = []
for p in PRIMES:
    k = 1
    while p ** (k + 1) <= B:
        k += 1
    SCHEDULE.extend([p] * k)


def measure(b):
    """True if base b FAILS on N (first nontrivial gcd event is N itself)."""
    acc = pow(b, 1, N)
    for step, p in enumerate(SCHEDULE):
        acc = pow(acc, p, N)
        # After the first multiplication (the first squaring) and every one after it,
        # check the gcd; the base's fate is decided by the FIRST nontrivial value.
        d = math.gcd(acc - 1, N)
        if d == N:
            return True
        if d > 1:
            return False
    return True  # gcd never left 1: no factor found, the base fails


def order(b, q):
    o, x, k = None, b % q, 1
    while x != 1:
        x = x * b % q
        k += 1
    return k


def key_corrected(o):
    if o in (1, 2):
        return (2, 1)
    fs = {}
    d = 2
    while d * d <= o:
        while o % d == 0:
            fs[d] = fs.get(d, 0) + 1
            o //= d
        d += 1
    if o > 1:
        fs[o] = fs.get(o, 0) + 1
    for l, e in fs.items():
        if l ** e > B:
            return None
    lmax = max(fs)
    return (lmax, fs[lmax])


def key_old(o):
    if o == 1:
        return (0, 0)
    return key_corrected(o)


def predict(b, key_fn):
    keys = [key_fn(order(b, q)) for q in FACTORS]
    exposed = [k for k in keys if k is not None]
    if not exposed:
        return True
    return all(k == min(exposed) for k in keys)


valid = [b for b in range(2, N - 1) if math.gcd(b, N) == 1]
meas = {b: measure(b) for b in valid}
pred_new = {b: predict(b, key_corrected) for b in valid}
pred_old = {b: predict(b, key_old) for b in valid}

n_meas = sum(meas.values())
n_new = sum(pred_new.values())
n_old = sum(pred_old.values())
dis_new = [b for b in valid if meas[b] != pred_new[b]]
dis_old = [b for b in valid if meas[b] != pred_old[b]]

print(f"valid bases: {len(valid)}")
print(f"measured failures: {n_meas} ({n_meas/len(valid):.4f})")
print(f"corrected model:   {n_new} ({n_new/len(valid):.4f}), disagreements: {len(dis_new)}")
print(f"old model:         {n_old} ({n_old/len(valid):.4f}), disagreements: {len(dis_old)} {dis_old}")
