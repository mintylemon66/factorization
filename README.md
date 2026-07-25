# factorization

Code for my paper and ICM 2026 poster, "Performance of Classical Factorization Algorithms for Structured Integers" (Isabella Li, Lexington High School; poster in Section 3, Number Theory, Philadelphia, July 2026).

The project compares classical factoring methods on three families of inputs: Mersenne numbers (2^i - 1 for i from 2 to 107), Fermat numbers (2^(2^n) + 1), and random 20 and 30 digit integers.

## Which Pollard method?

The Pollard timings in the original tables were produced by `3_pollards_rho.py`, which implements Pollard's rho. Early drafts and the conference abstract called that column Pollard's p-1, which is a different algorithm. I caught the mix-up while preparing the poster. The revised poster labels the column correctly and adds a separate column for an actual p-1 implementation, `4_pollards_p_minus_1.py` (B = 100,000, bases 2, 3 and 5, single stage, with checkpointed backtracking).

Re-running everything turned up something I did not expect. On the Mersenne inputs, the two Pollard methods fail on completely different numbers:

- rho fails exactly when i is a multiple of 20 (i = 20, 40, 60, 80, 100). Those inputs are divisible by 25, since the order of 2 mod 25 is 20, and rho with f(x) = x^2 + 1 cannot split 25: every starting value falls into the cycle 1 -> 2 -> 5, whose elements are distinct mod 5, so the gcd always comes back as 25 itself.
- p-1 fails exactly for i = 11, 29, 101. Every prime factor of 2^p - 1 is congruent to 1 mod 2p, so the factors tend to surface at the same multiplication and the gcd jumps from 1 straight to n.

Each method handles the other's failures instantly.

## Files

- `1_trial_division.py`: trial division (divide by 2, then odd candidates), runs to a full factorization
- `2_fermat.py`: Fermat's difference of squares method, using Decimal square roots
- `3_pollards_rho.py`: Pollard's rho (x^2 + 1, Floyd cycle detection, 8 restarts of 200,000 iterations) plus a recursive full factorization wrapper with a Miller-Rabin check. This file produced the original Pollard columns
- `4_pollards_p_minus_1.py`: Pollard's p-1, added July 2026. Produced the p-1 column on the revised poster
- `qs_fixed_minimal.py`: a minimal Quadratic Sieve. B is chosen adaptively from the input size (about 100 at 20 digits, about 400 at 30), one sieve interval of length about B^3, Tonelli-Shanks roots, Gaussian elimination over F2
- `20digit.txt`, `30digit.txt`: 500 random integers each, the datasets behind Tables 4 and 5 of the paper
- `list.py`: random integer dataset generator
- `randomsemiprime.py`: generates balanced 30 digit semiprimes (two random 15 digit primes), used for side experiments

## Timing methodology

Apple M3, CPU time via `time.process_time`, one run per input. Runtimes are heavy-tailed, so the paper reports medians and IQRs rather than means. Medians near 0.001 ms are at the timer's resolution floor, not real runtimes.

To reproduce a run: each numbered script has its timing loop at the bottom. Paste a dataset in as `my_list` (for example from `20digit.txt`) and run the script with python3. It prints one runtime in ms per input.

Two implementation notes for anyone reading closely: the QS gives up if a single sieve interval does not produce enough smooth relations (no retry, no multiple polynomials), and it proceeds at exactly |P| relations where the dependency guarantee needs |P| + 1. Both are deliberate simplifications of a baseline implementation.
