import math, json, importlib.util
BASE = './'
spec = importlib.util.spec_from_file_location('pm1', BASE + '4_pollards_p_minus_1.py')
pm1 = importlib.util.module_from_spec(spec); spec.loader.exec_module(pm1)

B_BOUND = 100000
def factorize_small(o):
    fs = {}; d = 2
    while d*d <= o:
        while o % d == 0: fs[d] = fs.get(d, 0) + 1; o //= d
        d += 1
    if o > 1: fs[o] = fs.get(o, 0) + 1
    return fs
def exposure_key(o):
    fs = factorize_small(o)
    if not fs: return (0, 0)
    for l, e in fs.items():
        if l**e > B_BOUND: return None
    lmax = max(fs)
    return (lmax, fs[lmax])
def order_small(b, q):
    o, x = 1, b % q
    while x != 1: x = x*b % q; o += 1
    return o

n = 2047
pred = meas = tot = 0
for b in range(2, n - 1):
    if math.gcd(b, n) != 1: continue
    tot += 1
    k1 = exposure_key(order_small(b, 23))
    k2 = exposure_key(order_small(b, 89))
    ex = [k for k in (k1, k2) if k is not None]
    p_fail = (not ex) or all(k == min(ex) for k in (k1, k2))
    if p_fail: pred += 1
    if pm1.p_minus_1(n, bases=(b,)) is None: meas += 1
print(json.dumps({'n_bases': tot,
                  'predicted_fail_fraction': round(pred/tot, 4),
                  'measured_fail_fraction': round(meas/tot, 4)}))
