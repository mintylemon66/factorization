import math, re, sys, time, json, signal
sys.setrecursionlimit(100000)
BASE = './'

# load her rho factor() wrapper
src = open(BASE + '3_pollards_rho.py').read()
ns = {}
exec(src[:src.index('my_list=')], ns)
factor_full, rho_once, is_pp = ns['factor'], ns['pollards_rho_once'], ns['is_probable_prime']

# load p-1 implementation
import importlib.util
spec = importlib.util.spec_from_file_location('pm1', BASE + '4_pollards_p_minus_1.py')
pm1mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(pm1mod)

class Timeout(Exception): pass
def _al(s, f): raise Timeout()
signal.signal(signal.SIGALRM, _al)

out = {}

# ---- A: which prime squares does x^2+1 split? exhaustive ----
def sieve(n):
    bs = [True]*(n+1); bs[0]=bs[1]=False
    for p in range(2, int(n**0.5)+1):
        if bs[p]: bs[p*p::p] = [False]*len(bs[p*p::p])
    return [i for i,b in enumerate(bs) if b]

def splits_square(p, c=1):
    m = p*p
    for s in range(m):
        x = y = s
        for _ in range(4*m):
            x = (x*x + c) % m
            y = (y*y + c) % m; y = (y*y + c) % m
            g = math.gcd(abs(x-y), m)
            if g == p: return True
            if x == y: break
    return False

unsplit = [p for p in sieve(100) if not splits_square(p)]
out['A_unsplittable_prime_squares_x2plus1'] = unsplit
print("A:", unsplit, flush=True)

# ---- B: which c in 1..24 split 25? ----
bad_c = [c for c in range(1, 25) if not splits_square(5, c)]
out['B_c_values_that_cannot_split_25'] = bad_c
print("B:", bad_c, flush=True)

# ---- C: sync prediction vs experiment ----
def full_factor(n, cap=600):
    signal.alarm(cap)
    try:
        fs = factor_full(n); signal.alarm(0); return sorted(fs)
    except Timeout:
        return None
    finally:
        signal.alarm(0)

KNOWN = {101: [7432339208719, 341117531003194129]}  # known factorization of 2^101 - 1

def mult_order(b, q):
    # order of b mod q via factoring q-1
    fs = full_factor(q-1, cap=300)
    if fs is None: return None
    o = q - 1
    for f in set(fs):
        while o % f == 0 and pow(b, o//f, q) == 1:
            o //= f
    return o

B_BOUND = 100000
def exposure_key(order):
    # her staged exponentiation: primes ascending, prime power l^e as e single multiplications.
    # factor q exposed when accumulated exponent divisible by its order; completing step is
    # the e-th multiplication of the LARGEST prime l in order (earlier primes done first).
    # unexposed if any prime power in order exceeds B (not powersmooth).
    fs = {}
    o = order
    d = 2
    while d*d <= o:
        while o % d == 0: fs[d] = fs.get(d,0)+1; o //= d
        d += 1
    if o > 1: fs[o] = fs.get(o,0)+1
    # same boundary correction as in stalls_density_repunits.py: the first gcd
    # check follows the first squaring, so orders 1 and 2 share a key
    if not fs: return (2, 1)
    for l, e in fs.items():
        if l**e > B_BOUND: return None
    lmax = max(fs)
    return (lmax, fs[lmax])

comp_exps = [p for p in sieve(107) if not is_pp((1<<p)-1) and p > 2]
predictions, factor_cache = {}, {}
for p in comp_exps:
    n = (1 << p) - 1
    fs = KNOWN.get(p) or full_factor(n)
    if fs is None:
        predictions[p] = 'factorization_timeout'; continue
    factor_cache[p] = fs
    fails_all = True
    for b in (2, 3, 5):
        keys = []
        inconclusive = False
        for q in sorted(set(fs)):
            o = mult_order(b, q)
            if o is None:
                inconclusive = True; break
            keys.append(exposure_key(o))
        if inconclusive:
            predictions[p] = 'inconclusive_order'; fails_all = None; break
        exposed = [k for k in keys if k is not None]
        if not exposed:
            base_fails = True
        else:
            kmin = min(exposed)
            base_fails = all(k == kmin for k in keys)  # every factor exposed at same first event -> gcd jumps to n
        if not base_fails:
            fails_all = False; break
    if fails_all is not None:
        predictions[p] = 'FAIL' if fails_all else 'success'
    print(f"C: p={p} predicted {predictions[p]}", flush=True)

out['C_predicted_failures'] = sorted([p for p, v in predictions.items() if v == 'FAIL'])
out['C_observed_failures'] = [11, 29, 101]
out['C_notes'] = {str(k): v for k, v in predictions.items() if v == 'factorization_timeout'}
print("C predicted:", out['C_predicted_failures'], "observed: [11, 29, 101]", flush=True)

# ---- C2: base-failure density on 2047, predicted vs measured ----
n = 2047
pred_fail = 0; tot = 0
for b in range(2, n-1):
    if math.gcd(b, n) != 1: continue
    tot += 1
    # compute orders directly (23, 89 tiny)
    def order_small(b, q):
        o, x = 1, b % q
        while x != 1: x = x*b % q; o += 1
        return o
    o1, o2 = order_small(b, 23), order_small(b, 89)
    k1, k2 = exposure_key(o1), exposure_key(o2)
    if k1 == k2: pred_fail += 1
out['C2_2047_predicted_fail_fraction'] = round(pred_fail/tot, 4)

meas_fail = 0; tot2 = 0
for b in range(2, n-1):
    if math.gcd(b, n) != 1: continue
    tot2 += 1
    if pm1mod.p_minus_1(n, bases=(b,)) is None: meas_fail += 1
out['C2_2047_measured_fail_fraction'] = round(meas_fail/tot2, 4)
print("C2:", out['C2_2047_predicted_fail_fraction'], "vs", out['C2_2047_measured_fail_fraction'], flush=True)

json.dump(out, open('blindspots_and_predictions_results.json','w'), indent=1)
print("DONE")
