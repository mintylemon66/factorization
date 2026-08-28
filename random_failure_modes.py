# Classify WHY p-1 fails on the random 20- and 30-digit datasets.
# Two possible failure modes:
#   SYNC  = gcd jumped from 1 straight to n (order synchronization), the
#           structured regime seen on 2^11-1, 2^29-1, 2^101-1
#   NOSMOOTH = no base ever produced any nontrivial gcd, i.e. no prime factor
#           has B-powersmooth p-1. Nothing to do with synchronization.
# Uses Isabella's implementation verbatim (4_pollards_p_minus_1.py), with the
# failure mode recorded.
import math, random, sys

B = 100_000
ICM = './'

def _sieve(limit):
    bs = [True] * (limit + 1)
    bs[0] = bs[1] = False
    for p in range(2, int(limit ** 0.5) + 1):
        if bs[p]:
            bs[p * p::p] = [False] * len(bs[p * p::p])
    return [i for i, b in enumerate(bs) if b]

_PRIMES = _sieve(B)
_PP = [(q, q ** int(math.log(B) / math.log(q))) for q in _PRIMES]

def is_probable_prime(a, k=10):
    if a < 2: return False
    for p in (2,3,5,7,11,13,17,19,23,29,31,37):
        if a % p == 0: return a == p
    d, r = a - 1, 0
    while d % 2 == 0: d //= 2; r += 1
    for _ in range(k):
        x = pow(random.randrange(2, a - 1), d, a)
        if x in (1, a - 1): continue
        for _ in range(r - 1):
            x = x * x % a
            if x == a - 1: break
        else: return False
    return True

def p_minus_1_diag(n, bases=(2,3,5), check_every=200):
    """Isabella's p_minus_1 plus a per-base record of how it ended."""
    modes = []
    if n % 2 == 0: return 2, ['even']
    for b in bases:
        g = math.gcd(b, n)
        if 1 < g < n: return g, modes + ['gcd_base']
        if g == n: modes.append('base_shares_n'); continue
        x = b
        ckpt, ckpt_i = x, 0
        overshoot = False
        for idx in range(0, len(_PP), check_every):
            for q, qe in _PP[idx:idx + check_every]:
                x = pow(x, qe, n)
            d = math.gcd(x - 1, n)
            if d == n:
                y = ckpt; stop = False
                for q, qe in _PP[ckpt_i:idx + check_every]:
                    e = 0
                    while q ** (e + 1) <= qe: e += 1
                    for _ in range(e + 1):
                        y = pow(y, q, n)
                        d2 = math.gcd(y - 1, n)
                        if 1 < d2 < n: return d2, modes + ['recovered_by_backtrack']
                        if d2 == n: stop = True; break
                    if stop: break
                overshoot = True
                break
            if d > 1: return d, modes + ['exposed']
            ckpt, ckpt_i = x, idx + check_every
        if overshoot:
            modes.append('SYNC')          # gcd hit n, backtracking could not separate
            continue
        d = math.gcd(x - 1, n)
        if 1 < d < n: return d, modes + ['exposed_final']
        modes.append('NOSMOOTH')          # never any nontrivial gcd for this base
    return None, modes

def load(path):
    txt = open(ICM + path).read()
    return [int(s) for s in txt[txt.index('[') + 1: txt.rindex(']')].split(',')]

for name in ('20digit.txt', '30digit.txt'):
    data = load(name)
    fails = []
    for n in data:
        if is_probable_prime(n):
            continue
        r, modes = p_minus_1_diag(n)
        if r is None:
            fails.append((n, modes))
    print(f'--- {name}: {len(data)} inputs, {len(fails)} p-1 failures')
    n_sync = sum(1 for _, m in fails if 'SYNC' in m)
    n_nosm = sum(1 for _, m in fails if 'SYNC' not in m)
    print(f'    synchronization failures (gcd 1 -> n): {n_sync}')
    print(f'    non-smooth failures (no factor exposed at all): {n_nosm}')
    for n, m in fails:
        print(f'    {n}  modes={m}')
    sys.stdout.flush()
