# Changelog

Accuracy figures are worst-case errors against mpmath on random exact inputs (300 points per region unless noted), double / x87 long double.

## std::format width for the whole quaternion; mixed-precision operators

* **Formatting (behavior change).**
  * Fill, align and width now apply to the whole `(a,b,c,d)` text, not to each component. This matches `operator<<` and the P2197 proposal for `std::complex`. The default alignment is right.
  * Sign, `#`, precision, `L` and type still apply to every component.
  * Nested `{}`/`{n}` width and precision are supported, and so is a multi-byte fill.
  * The `0` flag and non-floating types (`d`, `s`, `x`, …) are rejected with `std::format_error`.
  * Before: `std::format("{:>4}", q)` gave `(   1,   2,   3,   4)`. It now gives `(1,2,3,4)`, because a width of 4 is narrower than the text. Use `{:>20}` to pad the whole quaternion.
* **Mixed-precision operators.**
  * `+ - * /` between `quaternion<T>` and a floating-point `U ≠ T` or `std::complex<U>`, in either order, return `quaternion<common_type_t<T, U>>`, computed in that type. `==` against those types also compares in the common type.
  * Before: `quaternion<float> * 1.0` narrowed `1.0` to `float` silently. It is now a `quaternion<double>`.
  * Integers still convert to `T`.
  * `std::complex` has no mixed operators at all. The promotion follows the floating-point rule of [cmplx.over] (which applies to `pow` and friends). README explains the choice.
* **Fix (found by the manual fuzz run for this change, seed 2871196519):** `cyl_neumann`/`cyl_hankel_1` at large order and tiny |z| returned NaN when H⁽¹⁾_{ν+1} overflowed but the result did not. Example: ν ≈ 75, |z| ≈ 4.4e-3, Y ≈ −1.24e306. J, which is then below the smallest normal number, is now returned as 0, and Y/H⁽¹⁾ are now exact to < 1e-12 (before: NaN). Test `neumann_near_overflow`.
* **Fuzz script:**
  * Results within 2¹² of the largest finite value are skipped. Another Bessel intermediate still overflows there: Y_ν at ν = 84.5, |z| ≈ 0.014 gives −inf for Im Y ≈ −1.87e307 (seed 1304285165, clang-18/libc++). This is a documented limit.
  * Hurwitz ζ(s, a) with Re a < 0 is now measured against the largest term |a + n|^−s, because the terms cancel at its real zeros. One example from seed 1842008100: s = 0.83, a = −5.41, ζ ≈ 5.9e-4 measured 360 ε relative, which is about 1 ε of the largest term.
  * The Airy bound is now 128 ε, because the condition number grows like |z|^{3/2}. The largest error seen was 69.6 ε in long double at z ≈ −18.6 + 5.8i.
* **Tests:**
  * `std_format`: 23 new checks across 3 types, including error cases.
  * `mixed_precision`: result types via `static_assert`, values, complex on either side, equality, and `constexpr` use.

## Nightly mpmath fuzz job; fixes for ζ near 0, Li_s near integer s, and W₀

* **Fuzz job.** `.github/workflows/fuzz.yml` runs nightly (cron 06:17 UTC) and on `workflow_dispatch` (optional `seed`, `samples`). It is a separate, non-gating workflow, with GCC 14/libstdc++ and Clang 18/libc++.
  * `fuzz/fuzz.py` samples 22 special functions over several regions: real-axis points, |z| from 10⁻¹² to 400, Bessel orders from −5.5 to 100, ζ near 0 and 1, and Li_s near integer s. Points go on random ±i/±j/±k axes.
  * The evaluator `fuzz/fuzz_eval.cpp` computes the values, and the script compares them with mpmath against per-function bounds. The seed is logged, and the report goes to the job summary and to the `fuzz-report-*` artifact.
  * `quaternion_fuzz_eval` is built with the tests, so the evaluator stays compiling on every CI toolchain.
* **Bugs it found** (before → after; worst error in ε over 400 fuzz samples, seed 11, double / long double):
  * **`riemann_zeta` next to s = 0** (|s| ≲ 10⁻³): the functional equation used ζ(1 − s) with 1 − s rounded, so the error was ε/|s|. Examples: 2.6e6 ulp at s = −10⁻⁷, and 4e11 ulp at s = −10⁻¹². It now uses the Taylor series about 0 for |s| < 0.25. Fuzz: 1.8e11 / 1.4e11 → 11.8 / 9.3.
  * **`polylog(s, z)` for non-integer s next to a positive integer n:**
    * The k = n−1 term ζ(1+δ)μⁿ⁻¹/(n−1)! and Γ(1−s)(−μ)^{s−1} both have 1/δ poles. They are now combined analytically (Stieltjes constants and a series for log(πδ/sin πδ) + log Γ(n) − log Γ(n+δ)).
    * The log-series stopped early at a term next to a trivial zero of ζ(s−k). ζ there now comes from the reflection formula with the exact offset from −2m.
    * Fuzz: 3.5e13 / 7.1e16 → 100 / 109. At s = 3.00000001, z = −0.8: 6.9e13 / 1.2e15 → < 64.
  * **`lambert_w` (branch 0) for 1 < |z| ≤ 3:** the start log z − log log z could send Halley's iteration to W₋₁ or W₁ (e.g. z = 0.9668 + 0.2678i gave −1.59 + 4.65i). On a grid over 1 < |z| ≤ 6, 407 of 27,692 points were affected. It now starts from log(1+z), or from the branch-point series near z = −1, and checks that the result lies in the range of W₀, retrying from other starts if it does not. Fuzz: 3.6e16 / 7.5e19 → 1.0 / 1.1. Grid: 407 wrong → 0 of 32,800.
* **Not library bugs:** the fuzz also flagged three other things:
  * sph_bessel/sph_neumann on the negative real axis;
  * K/J results that under- or overflow `double`;
  * H⁽¹⁾ for Im z ≳ 40, where mpmath's `hankel1` cancels in J + iY and loses digits, so the script now raises the working precision by the digits lost.
  The first case was a reference issue: mpmath's √(π/2z)·J_{n+½} takes the cut, while j_n and y_n are single-valued. The script now uses the parity relations there and skips results outside the normal range of the type.
* **Tests:** `fuzz_regressions_zeta_polylog_lambert` (17 mpmath references).

## Faster Bessel functions at small |z|; relative accuracy next to the zeros of J

* **J_ν, I_ν, j_n and integer-order J_n at |z| ≤ 16** now use the ascending series. The sum S = Σ(∓z²/4)ᵏ/(k!(ν+1)ₖ) is carried in double-length arithmetic, which absorbs its cancellation of up to e^|z|. The prefactor (z/2)^ν/Γ(ν+1) is formed as a double-length integer power times (z/2)^frac(ν). For large order, it is built as exp(ν log(z/2) − log Γ(ν+1)) with a double-length exponent. Negative non-integer ν uses the series directly, and negative integers use J₋ₙ = (−1)ⁿJₙ and I₋ₙ = Iₙ. These functions used to go through K (the double-exponential Laplace integral) and the J/I continued fraction, or, for integer n, Miller's recurrence, which is only absolutely accurate.
* **K_ν at small and moderate |w|.** Temme's series is used for |w| ≤ 2, with 1/Γ(1 ± a) and Temme's γ₁, γ₂ taken from a 33-term Taylor series of 1/Γ(1+x). Steed's CF2 is used for 2 < |w| < 8, and wherever Hankel's expansion does not converge. Half-integer orders use closed forms. All of these start from the base order a = ν − round(ν), |a| ≤ ½, and recur upward. Both sums are compensated. The Laplace integral is kept only as a fallback if CF2 does not converge. Y, H⁽¹⁾, H⁽²⁾, K, the Airy functions and sph_neumann all gain from this.
* **x87 `long double`:** `two_prod` uses Dekker's exact product (Veltkamp split) instead of `std::fma`, which is a slow software routine on x87. This speeds up every double-length kernel in `long double`.
* **Tests:** `bessel_near_zeros_and_small_k` checks 92 mpmath references for J at, and next to, the zeros j_ν,k ≤ 16 (ν = 0, ¼, 1, 2.5, 3, 7). The integer-order kernel is checked to agree with them. It also checks 36 references for K at |w| ≤ 8.
* **Benchmarks:** new cases for J, I, j_n, Y and K at |z| ≈ 0.5–10.

Accuracy, before → after. Errors are against mpmath at 811 random and zero-adjacent exact points with |z| ≤ 16, double / long double:

| | double | long double |
|---|---|---|
| J_ν, ν ∈ {0, 0.3, 1, …, 12, −0.3, −1.5, −2.5} | 13.9 → 4.7 ulp | 9.8 → 3.4 ulp |
| J_n (integer kernel) | 45 → 2.3 ulp | 35 → 4.2 ulp |
| j_n (sph_bessel) | 51 → 4.1 ulp | 11.6 → 3.9 ulp |
| I_ν | 7.7 → 4.4 ulp | 9.7 → 4.3 ulp |
| J at 10⁻⁸ or 10⁻¹² off a real zero | 1.6e13 → 14 ulp | 1.6e13 → 3.6 ulp |
| J at the double nearest a zero, absolute error | 3.6e-16 → 2.9e-28 | 1.7e-19 → 1.4e-34 |
| K_ν, \|w\| ≤ 8 (300 extra points, many near the imaginary axis) | 18.4 → 16.6 ulp | 18.7 → 17.5 ulp |
| Y_ν | 10.8 → 11.2 ulp | 12.1 → 18.3 ulp |

At 8 ≤ |z| ≤ 60, K and Y are unchanged (35 / 33 ulp and 26 / 37 ulp).

Timings, before → after, in ns per call. Measured with GCC 14 -O2 on x86-64 with `bench/bench.cpp`; quaternion arguments:

| case | double | long double |
|---|---|---|
| J_2.5, \|z\| ≈ 0.5 | 24,437 → 806 (30×) | 143,899 → 3,491 (41×) |
| J_2.5, \|z\| ≈ 2 | 24,290 → 1,162 (21×) | 141,365 → 4,803 (29×) |
| J_2.5, \|z\| ≈ 6 | 16,616 → 1,557 (11×) | 94,472 → 6,754 (14×) |
| J_0.3, \|z\| ≈ 10 | 23,783 → 1,765 (13×) | 149,831 → 8,841 (17×) |
| I_1.5, \|z\| ≈ 2 | 26,992 → 1,091 (25×) | 152,893 → 4,708 (32×) |
| I_0.3, \|z\| ≈ 6 | 28,906 → 1,984 (15×) | 146,495 → 6,530 (22×) |
| sph_bessel(2), \|z\| ≈ 2 | 25,143 → 1,321 (19×) | 144,530 → 4,834 (30×) |
| Y_1.5, \|z\| ≈ 2 | 26,329 → 644 (41×) | 147,829 → 2,610 (57×) |
| Y_0.3, \|z\| ≈ 6 | 31,287 → 2,179 (14×) | 147,050 → 10,085 (15×) |
| K_0.3, \|z\| ≈ 2 | 27,864 → 3,015 (9×) | 150,571 → 14,602 (10×) |
| K_0.5, \|z\| ≈ 2 | 24,671 → 256 (96×) | 145,161 → 1,272 (114×) |
| J_3 (integer), \|z\| ≈ 6 | 893 → 1,485 (0.6×; relative accuracy near zeros) | 4,926 → 6,476 (0.8×) |
| large \|z\| / large ν cases | unchanged (±10%) | unchanged |

## Real-axis special-function coverage on libc++; Clang 19 in CI

* `tests/realaxis_refs.inc` holds 366 mpmath references (40 digits), generated by `tests/gen/gen_realaxis.py`. It uses the `std::` conventions: no Condon–Shortley phase in `assoc_legendre`, modulus `k` for the elliptic integrals, and (1 − ν sin²θ) in Π. The inputs are those of the `std::` special-math comparisons, which only run on libstdc++. Inexact decimals are replaced by nearby binary fractions.
  * Coverage: J, Y, I, K of real order; integer-order J and I; spherical Bessel; Hermite, Laguerre, associated Laguerre, Legendre, associated and spherical Legendre; Ei; the incomplete and complete elliptic integrals of all three kinds; ζ; and β.
  * The test `real_axis_mpmath_all_libraries` checks them on every standard library, at 64 ε. Polynomials use absolute error for |value| < 1, and J/Y for x > 1 use error relative to the envelope √(2/(πx)).
  * The worst errors are 22 ε (double) and 32 ε (long double); float results are exact after rounding.
* **Bug fix found by these tables.** Since the previous release, `cyl_bessel_i(ν, x)` for real x ≥ 16 (on the axis, with ν small enough for the asymptotic expansion) returned a wrong value (e.g. `I_0(40)`). `sqrt(π/(2w))` flipped the sign of the zero imaginary part of −w, which took the wrong branch. It is now `sqrt(π/2)/sqrt(w)`.
* **CI.** Clang 19 (from apt.llvm.org) on libstdc++ and libc++, Debug and Release, added to the matrix.

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
