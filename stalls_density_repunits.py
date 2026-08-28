import math, sys, time, json, signal, importlib.util
BASE = './'
src = open(BASE + '3_pollards_rho.py').read()
ns = {}
exec(src[:src.index('my_list=')], ns)
factor_full, is_pp = ns['factor'], ns['is_probable_prime']
spec = importlib.util.spec_from_file_location('pm1', BASE + '4_pollards_p_minus_1.py')
pm1 = importlib.util.module_from_spec(spec); spec.loader.exec_module(pm1)

class Timeout(Exception): pass
def _al(s, f): raise Timeout()
signal.signal(signal.SIGALRM, _al)
def with_cap(fn, cap, *a):
    signal.alarm(cap)
    try: r = fn(*a); return r
    except Timeout: return None
    finally: signal.alarm(0)

out = {}

# ---- E5: does x^2+2 repair rho on the five failed inputs? ----
import random
def rho_c_once(n, c):
    if n % 2 == 0: return 2
    for _ in range(8):
        x = y = random.randrange(2, n - 1); d = 1
        for _ in range(200000):
            x = (x*x + c) % n
            y = (y*y + c) % n; y = (y*y + c) % n
            d = math.gcd(abs(x - y), n)
            if 1 < d < n: return d
            if d == n: break
    return None
stall_box = []
def factor_fixed_c(n, c):
    if n == 1: return []
    if is_pp(n): return [n]
    d = None
    tries = 0
    while d is None:
        d = rho_c_once(n, c); tries += 1
        if tries > 6:
            stall_box.append(n)
            raise Timeout()
    return factor_fixed_c(d, c) + factor_fixed_c(n // d, c)
def factor_rand_c(n):
    if n == 1: return []
    if is_pp(n): return [n]
    d = None
    while d is None: d = rho_c_once(n, random.randrange(1, n - 1))
    return factor_rand_c(d) + factor_rand_c(n // d)

e5 = {}
for i in (20, 40, 60, 80, 100):
    n = (1 << i) - 1
    stall_box.clear()
    fs1 = with_cap(factor_fixed_c, 45, n, 2)
    fs2 = with_cap(factor_rand_c, 90, n)
    ok2 = fs2 is not None and math.prod(fs2) == n and all(is_pp(f) for f in fs2)
    e5[i] = {'fixed_c2': 'success' if fs1 is not None else ('stalls at ' + str(stall_box[0]) if stall_box else 'stalls (timeout)'),
             'randomized_c': 'success' if ok2 else 'fails'}
    print("E5:", i, e5[i], flush=True)
out['E5_c_variants_on_rho_failures'] = e5
# verify the 15 phase-lock claim directly
out['E5_rho_c2_on_15'] = 'splits' if with_cap(rho_c_once, 20, 15, 2) not in (None,) else 'cannot split (phase-locked)'
print("15 under c=2:", out['E5_rho_c2_on_15'], flush=True)

# ---- shared machinery: exposure keys ----
B_BOUND = 100000
def factorize_small(o):
    fs = {}; d = 2
    while d*d <= o:
        while o % d == 0: fs[d] = fs.get(d, 0) + 1; o //= d
        d += 1
    if o > 1: fs[o] = fs.get(o, 0) + 1
    return fs
def exposure_key_from_order(o):
    fs = factorize_small(o)
    # Corrected 2026-08-17: the first gcd check happens after the first squaring,
    # not at the base itself, so orders 1 and 2 are first detected at the same
    # check and must share a key. Giving order 1 a class of its own (the old
    # (0,0)) mispredicts exactly two of the 1934 valid bases on 2047.
    if not fs: return (2, 1)
    for l, e in fs.items():
        if l**e > B_BOUND: return None
    lmax = max(fs)
    return (lmax, fs[lmax])
def order_mod(b, q, qm1_factors):
    o = q - 1
    for f in set(qm1_factors):
        while o % f == 0 and pow(b, o // f, q) == 1: o //= f
    return o

# ---- C3: predicted vs measured base-failure density, M29 and 2047 spot-check family ----
def density(n, factors, sample_bases):
    qm1f = {}
    okq = True
    for q in factors:
        fs = with_cap(factor_full, 300, q - 1)
        if fs is None: okq = False; break
        qm1f[q] = fs
    if not okq: return None
    pred_fail = tot = 0
    for b in sample_bases:
        if math.gcd(b, n) != 1: continue
        tot += 1
        keys = []
        bad = False
        for q in factors:
            k = exposure_key_from_order(order_mod(b, q, qm1f[q]))
            keys.append(k)
        exposed = [k for k in keys if k is not None]
        if not exposed: pred_fail += 1; continue
        kmin = min(exposed)
        if all(k == kmin for k in keys): pred_fail += 1
    meas_fail = tot2 = 0
    for b in sample_bases:
        if math.gcd(b, n) != 1: continue
        tot2 += 1
        if pm1.p_minus_1(n, bases=(b,)) is None: meas_fail += 1
    return {'predicted': round(pred_fail / tot, 4), 'measured': round(meas_fail / tot2, 4), 'n_bases': tot}

m29 = (1 << 29) - 1
fs29 = with_cap(factor_full, 120, m29)
if fs29:
    out['C3_M29_density'] = density(m29, sorted(set(fs29)), range(2, 1502))
    print("C3 M29:", out['C3_M29_density'], flush=True)

# ---- R: repunits, out-of-sample prediction of the p-1 failure set ----
rep = {}
for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41):
    R = (10**p - 1) // 9
    if is_pp(R):
        rep[p] = {'status': 'prime, skipped'}; print("R:", p, "prime", flush=True); continue
    fs = with_cap(factor_full, 300, R)
    if fs is None:
        rep[p] = {'status': 'factorization timeout'}; print("R:", p, "timeout", flush=True); continue
    qs = sorted(set(fs))
    qm1f, ok = {}, True
    for q in qs:
        f2 = with_cap(factor_full, 300, q - 1)
        if f2 is None: ok = False; break
        qm1f[q] = f2
    if not ok:
        rep[p] = {'status': 'order factorization timeout'}; continue
    fails_all = True
    for b in (2, 3, 5):
        keys = [exposure_key_from_order(order_mod(b, q, qm1f[q])) for q in qs]
        exposed = [k for k in keys if k is not None]
        if not exposed: base_fails = True
        else:
            kmin = min(exposed)
            base_fails = all(k == kmin for k in keys)
        if not base_fails: fails_all = False; break
    predicted = 'FAIL' if fails_all else 'success'
    r = with_cap(pm1.p_minus_1, 120, R)
    observed = 'FAIL' if r is None else 'success'
    rep[p] = {'predicted': predicted, 'observed': observed, 'match': predicted == observed}
    print("R:", p, rep[p], flush=True)
out['R_repunit_out_of_sample'] = rep

json.dump(out, open('stalls_density_repunits_results.json', 'w'), indent=1)
print("DONE")
