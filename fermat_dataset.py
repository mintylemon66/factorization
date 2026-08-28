# Times every method on the Fermat numbers F_0 through F_7, one run per input
# with time.process_time, like the other timing loops in this repo. This is the
# script behind the Fermat-number tables in the paper.
#
# Not attempted (would run for years or worse): trial division on F_7, Fermat's
# method on F_6, rho on F_7 (about 47 CPU-days expected under the 8 x 200,000
# restart budget), QS on F_7 (the sieve interval alone exceeds memory).
# Runs that finish without a factorization print f: p-1 on F_7 (no factor's
# order is 10^5-powersmooth), QS on F_5 and F_6 (too few smooth relations),
# Fermat's method on F_7 (decimal precision error at the first step).
#
# Run from the repo root: python3 fermat_dataset.py

import ast
import time


def load(fname, names):
    """Load functions from a repo script without running its timing driver."""
    tree = ast.parse(open(fname).read())
    keep = [n for n in tree.body
            if isinstance(n, (ast.Import, ast.ImportFrom, ast.FunctionDef, ast.ClassDef))
            or (isinstance(n, ast.Assign) and len(ast.dump(n)) < 2000)]
    ns = {}
    exec(compile(ast.Module(body=keep, type_ignores=[]), fname, 'exec'), ns)
    return [ns[x] for x in names]


(trial_division,) = load('1_trial_division.py', ['trial_division'])
(fermat_factor,) = load('2_fermat.py', ['fermat_factor'])
(rho_factor,) = load('3_pollards_rho.py', ['factor'])
p_minus_1, is_probable_prime = load('4_pollards_p_minus_1.py', ['p_minus_1', 'is_probable_prime'])
(QS,) = load('qs_fixed_minimal.py', ['QS'])

F = [2 ** (2 ** n) + 1 for n in range(8)]


def timed(fn, x):
    start = time.process_time()
    result = fn(x)
    ms = (time.process_time() - start) * 1000
    return ms, result


PLAN = [
    ('trial division', trial_division, range(7)),
    ("Fermat's method", fermat_factor, range(6)),
    ('rho (factor wrapper)', rho_factor, range(7)),
    ('p-1 (with primality pre-check)', lambda x: x if is_probable_prime(x) else p_minus_1(x), range(8)),
    ('QS', QS, range(7)),
]

for name, fn, idx in PLAN:
    print(f"== {name} ==")
    for n in idx:
        try:
            ms, r = timed(fn, F[n])
            tag = 'f' if r is None else str(r)[:60]
            print(f"  F_{n}: {ms:.3f} ms -> {tag}")
        except Exception as e:
            print(f"  F_{n}: f ({type(e).__name__}: {e})")
