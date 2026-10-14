# Changelog

Accuracy figures are worst-case errors against mpmath on random exact inputs (300 points per region unless noted), double / x87 long double.

## Benchmarks; faster Bessel (large |z|, large ν) and ζ (large |Im s|)

* **Benchmark suite.** `bench/bench.cpp` (CMake option `QUATERNION_BUILD_BENCHMARKS`, default OFF) is self-contained, with no dependencies. CI runs it in a non-gating job and uploads `benchmark-results` (`bench.json`, `bench.txt`).
* **Bessel functions.**
  * Hankel's asymptotic expansion is used for K, H⁽¹⁾, H⁽²⁾, J, Y and I at |z| ≥ 8–16, wherever its terms reach ε/4 while decreasing.
  * For ν ≥ 2, K is recurred upward from frac(ν) (Temme).
  * Integer-order J uses the expansion for |z| ≥ 16 when n² < |z|.
* **ζ.**
  * The Euler–Maclaurin sum uses N = 0.3|s| + 30 (was 0.6|s| + 20), with up to 80 Bernoulli corrections (B₂ⱼ/(2j)! via ζ(2j) for j > 20).
  * log k comes from a prime sieve in double-length arithmetic, so only primes need a logarithm.
  * Riemann–Siegel was not used: with the few correction terms that are practical for general complex s, it cannot reach ulp accuracy at |Im s| ≲ 10⁶.

Timings, before → after, in ns per call. Measured with GCC 14 -O2 on x86-64, quaternion arguments on a rotated axis:

| case | double | long double |
|---|---|---|
| J_2.5, \|z\| ≈ 30 | 25,071 → 404 (62×) | 141,714 → 2,007 (71×) |
| J_2.5, \|z\| ≈ 3000 | 38,161 → 421 (91×) | 208,242 → 2,710 (77×) |
| J_3 (integer), \|z\| ≈ 3000 | 46,636 → 594 (79×) | 282,786 → 4,016 (70×) |
| Y_1.5, \|z\| ≈ 300 | 30,541 → 422 (72×) | 162,379 → 2,552 (64×) |
| K_0.5, \|z\| ≈ 300 | 17,147 → 231 (74×) | 97,534 → 1,485 (66×) |
| I_1.5, \|z\| ≈ 300 | 33,026 → 243 (136×) | 164,481 → 1,398 (118×) |
| J_50, \|z\| ≈ 30 | 39,134 → 1,489 (26×) | 202,294 → 7,135 (28×) |
| J_200, \|z\| ≈ 300 | 248,090 → 2,959 (84×) | 1,107,100 → 12,283 (90×) |
| K_100, \|z\| ≈ 30 | 58,099 → 1,097 (53×) | 336,964 → 4,683 (72×) |
| J_2.5, \|z\| ≈ 2 (unchanged path) | 25,761 → 25,063 | 147,984 → 140,533 |
| ζ(½ + i·10²) | 151,127 → 43,302 (3.5×) | 3.40 ms → 0.99 ms (3.4×) |
| ζ(½ + i·10³) | 1.23 ms → 0.156 ms (7.9×) | 27.4 ms → 3.5 ms (7.8×) |
| ζ(½ + i·10⁴) | 11.8 ms → 1.01 ms (11.7×) | 215 ms → 23.1 ms (9.3×) |
| ζ(½ + i·10⁵) | 94.6 ms → 8.06 ms (11.7×) | 2108 ms → 138 ms (15×) |

Accuracy against mpmath (worst / p99 / median ulp, before → after, 100–150 random exact inputs per row; double, then long double). The Bessel rows all improve; ζ is unchanged within sampling noise (its long double worst case moves from 10.4 to 11 ulp):

| region | double | long double |
|---|---|---|
| J_ν, \|ν\| ≤ 3, 16 ≤ \|z\| ≤ 3000 | 720/411/5.8 → 7.1/4.5/1.7 | 4920/1170/6.2 → 18/5.6/1.8 |
| J_n, n ≤ 10, 16 ≤ \|z\| ≤ 3000 | 486/393/1.8 → 56/7.5/1.4 | 22300/561/1.7 → 38/8.3/1.4 |
| Y_ν, \|ν\| ≤ 3, large \|z\| | 4940/1990/5.8 → 15/9.3/1.8 | 1350/859/6.7 → 7.7/4.7/1.6 |
| K_ν, ν ≤ 3, large \|z\| | 561/233/2.8 → 6.3/5.8/1.5 | 195/60/2.8 → 5.9/5.5/1.8 |
| I_ν, ν ≤ 3, large \|z\| | 779/757/5.6 → 5.4/3.8/1.5 | 1300/126/4.7 → 16/8.7/1.6 |
| J_ν, 20 ≤ ν ≤ 250, \|z\| ≤ 400 | 1.2e15/4.3e8/376 → 649/272/44 | 4.4e9/5.4e8/419 → 242/201/61 |
| Y_ν, same | 6.5e15/8.9e12/346 → 294/239/41 | 2.4e15/3.6e9/324 → 325/227/35 |
| K_ν, same | 1.7e12/6.5e4/281 → 190/182/27 | 2.1e7/7.2e5/251 → 272/253/47 |
| I_ν, same | 7.5e9/6.4e3/318 → 433/219/33 | 8.8e12/3.2e6/347 → 524/229/50 |
| J_ν, K_ν, \|z\| ≤ 16 | unchanged | unchanged |
| ζ, 2 ≤ \|Im s\| ≤ 10⁵, −2 ≤ Re s ≤ 3 | 16.6/4.3/0.87 → 14.2/4.3/0.91 | 10.4/3.4/0.65 → 11/3.4/0.65 |

## C Annex G infinities, NaNs and signed zeros

The elementary functions, the arithmetic operators, `abs`, `arg`, `proj` and `polar` now follow C11/C23 Annex G. The special values are mapped to quaternions by axis reduction (see the README section "Infinities, NaNs and signed zeros").

* Every Annex G table entry, from `cexp` through `cacos` plus the derived `csin`/`ccos`/`ctan`/`casin`/`catan`, is reproduced on the `i` axis and lifted to arbitrary axes. This is tested for float, double and long double (`tests/test_annexg.cpp`; 1382 table points per type, each on six axes or with both signs of the `j`/`k` zeros).
* The library's own kernels are used, so the results are identical on libstdc++ and libc++. glibc's `cacosh(±0 + iNaN) = NaN + iπ/2` deviates from the table, and the library returns `NaN + iNaN`.
* `*` and `/` recover infinities (C11 G.5.1 `_Cmultd`/`_Cdivd`, generalized to four components). `/` is now `p·conj(q)/‖q‖²` with power-of-two scaling of `q`, so it no longer overflows for huge `q`.
* A −0 `i` component of a real argument now selects the lower side of a branch cut, as `std::complex` does: `log(−1 − 0i) = −πi`, `sqrt(−4 − 0i) = −2i`, `asin(2 − 0i) = π/2 − 1.3169…i`. Before, it was ignored.
* A zero component times an infinite factor gives a signed zero, not NaN. This affects every lifted function, including the special functions at their poles (e.g. `atan(±j)`).
* Finite-input fixes found by the sweep:
  * `sqrt` no longer cancels for `a < 0` with a tiny `|v|`.
  * `tan` and `tanh` use Kahan's complex tanh, so they no longer overflow to NaN for large `|a|` or `‖v‖`.
* A `-ffast-math` caveat is documented.

Finite inputs: the worst error against mpmath, as before → after, over 1000 random exact inputs per function. A fifth of the inputs are on the `i` axis, a fifth are real, and the magnitudes are log-uniform in 1e-3–1e6 (exp/sin/cos/sinh/cosh/pow: 1e-3–300). Figures are float / double / long double in ulp.

| function | before | after |
|---|---|---|
| sqrt | 4.1e6 / 3.2e14 / 1.9e16 | 2.4 / 2.8 / 2.7 |
| tan | NaN / NaN / NaN | 3.5 / 2.7 / 3.6 |
| tanh, \|a\| or ‖v‖ large | NaN / NaN / NaN | ≤ 3.6 (except ill-conditioned tan(‖v‖) for ‖v‖ ≳ 1e5) |
| tanh, median | 1.7e7 / 1.9e4 / 2.0 | ≈ 0 / ≈ 0 / ≈ 0 |
| asinh | 1.9 / 2.9 / 4.6 | 1.9 / 2.0 / 2.4 |
| div | 2.2 / 2.3 / 2.2 | 2.2 / 2.0 / 2.0 |
| exp log log10 sin cos sinh cosh asin acos atan acosh atanh mul pow | — | identical worst case (median within ±0.05 ulp) |

## Riemann zeta at large |Im s| (Euler–Maclaurin)

ζ(s) for |Im s| > 1 now uses Euler–Maclaurin summation in double-length arithmetic, with no cap on |Im s|. It was Borwein's η algorithm, capped at 6000 terms (|Im s| ≲ 6,700).

| region | before | after |
|---|---|---|
| −2 ≤ Re s ≤ 3, 1 ≤ Im s ≤ 1500 (200 pts) | 90 / 41 ulp | 19 / 14 ulp |
| Im s 1,500–5,000 (30 pts) | 240 / 64 ulp | 2.9 / 3.7 ulp |
| Im s 5,000–30,000 (30 pts; 10 pts for before) | 4e15 / 9e18 ulp (wrong) | 3.6 / 3.1 ulp |
| Im s 30,000–120,000 (30 pts) | wrong | 3.3 / 4.0 ulp |
| \|Im s\| ≤ 3, −30 ≤ Re s ≤ 30 (200 pts) | 6.3 ulp | 5.1 ulp |

## 6368066 — incomplete gamma Γ(a, z)

Backward-evaluated continued fraction, double-length z^a e^{−z} prefactor, and an adaptive choice between the series forms.

| region | before | after |
|---|---|---|
| ½ < a < 1 near the negative real axis | 150 / 169 ulp | 8.7 / 6.6 ulp |
| arg z ≈ 0.85π–0.95π | 103 / 130 ulp | 4.3 / 4.3 ulp |
| full plane, −8 ≤ a ≤ 8 | 80 / 91 ulp | 28 / 11 ulp |
| within 0.05π of the negative real axis, \|z\| ≤ 100 | 52 / 56 ulp | 18.5 / 21 ulp |
| a at or near 0, −1, −2, −4, −7 | 87 / 80 ulp | 8.1 / 13 ulp |
| 0.5 ≤ a ≤ 20, all arg z | 370 / 162 ulp | 43 / 15–20 ulp (within input conditioning) |
