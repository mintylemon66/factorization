# E5 recheck (2026-08-16): record the ACTUAL cofactor where fixed-c=2 rho stalls.
# An earlier version of stalls_density_repunits.py hard-coded the label "stalls (15 blind spot)"; that is
# wrong for i=60, where 9 | 2^60-1 and x^2+2 cannot split 9, so the stall is at 45.
# This script redoes only experiment E5, records the true stall cofactors, and
# updates the E5 entries of stalls_density_repunits_results.json (a .bak copy is saved first).
import json, math, random, shutil, os

HERE = os.path.dirname(os.path.abspath(__file__))

def is_probable_prime(n):
    if n < 2: return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0: return n == p
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True

def rho_c_once(n, c, iters=200000, restarts=8):
    if n % 2 == 0: return 2
    for _ in range(restarts):
        x = y = random.randrange(2, n - 1)
        for _ in range(iters):
            x = (x * x + c) % n
            y = (y * y + c) % n; y = (y * y + c) % n
            d = math.gcd(abs(x - y), n)
            if 1 < d < n: return d
            if d == n: break
    return None

def factor_fixed_c(n, c, tries_cap=6):
    # returns (primes, stalled_composites): every leaf is either a prime factor
    # or a composite that resisted tries_cap full rho attempts with this fixed c
    primes, stalls, stack = [], [], [n]
    while stack:
        m = stack.pop()
        if m == 1: continue
        if is_probable_prime(m):
            primes.append(m); continue
        for _ in range(tries_cap):
            d = rho_c_once(m, c)
            if d is not None: break
        else:
            stalls.append(m); continue
        stack += [d, m // d]
    return primes, stalls

if __name__ == '__main__':
    random.seed(2026)
    results = {}
    for i in (20, 40, 60, 80, 100):
        n = (1 << i) - 1
        primes, stalls = factor_fixed_c(n, 2)
        assert stalls, f'i={i}: unexpectedly factored completely'
        assert len(stalls) == 1, f'i={i}: multiple stalled composites {stalls}'
        assert math.prod(primes) * math.prod(stalls) == n
        results[str(i)] = stalls[0]
        print(f'i={i}: fixed c=2 stalls at {stalls[0]}', flush=True)

    path = os.path.join(HERE, 'stalls_density_repunits_results.json')
    shutil.copyfile(path, path + '.bak')
    with open(path) as f: data = json.load(f)
    for i, stall in results.items():
        data['E5_c_variants_on_rho_failures'][i]['fixed_c2'] = f'stalls at {stall}'
    data['E5_rho_c2_on_45'] = 'cannot split (mod-9 collapse onto {2,6} plus lockstep with mod-5)'
    data['E5_note'] = ('corrected 2026-08-16: original label hard-coded "15 blind spot" for all five; '
                       'i=60 actually stalls at 45 because 9 | 2^60-1 and x^2+2 never yields gcd 3 mod 9')
    with open(path, 'w') as f: json.dump(data, f, indent=1)
    print('stalls_density_repunits_results.json updated (backup saved as .bak)')
