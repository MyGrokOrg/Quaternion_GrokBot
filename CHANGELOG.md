# Changelog

Accuracy figures are worst-case errors against mpmath on random exact inputs (300 points per region unless noted), double / x87 long double.

## MSVC /fp:contract: the headers restore the includer's fp_contract state

* **Fix** (MSVC with `/fp:contract` or `/fp:fast` only; no change elsewhere, and no change to the library's results). The headers ended with `#pragma fp_contract(on)`, so a user's `#pragma fp_contract(off)` before the include was overridden for the rest of the translation unit. They now save the state with `#pragma float_control(push)` before their `fp_contract(off)` and restore it with `#pragma float_control(pop)` at their end. MSVC's `float_control` stack includes the contraction state: verified in CI on x64 (`/arch:AVX2 /fp:contract`) and ARM64 (`/fp:contract`).
  * Before → after, `test_fpcontract`, `a * b + c` defined after the include with `fp_contract(off)` before it: fused → not fused, on both x64 and ARM64. With the command-line state or `fp_contract(on)` before the include it is fused, before and after the include, in both versions.
  * The library's own code is still compiled without contraction: the Windows bit-identity jobs still find MSVC `/fp:contract` bit-identical to plain MSVC.
* **Tests.** New `test_fpcontract` (`tests/test_fpcontract.cpp`, `tests/fpcontract/off.cpp`, `tests/fpcontract/on.cpp`): three translation units (state from the command line, `fp_contract(off)` and `fp_contract(on)` before the include) compare `a * b + c` defined before and after the include, bit for bit, at a = 1 + 2⁻³⁰, b = 1 − 2⁻³⁰, c = −1 (0 unfused, −2⁻⁶⁰ fused). Under MSVC it also checks that the states take effect (`fp_contract(off)` not fused; with `/fp:contract` and an fma instruction the others fused). It runs on every toolchain; there are now 10 test suites.
* **CI.** New *windows* job *msvc ARM64 /fp:contract (Release)*; the *windows* jobs print `test_fpcontract`'s output. CI has 25 jobs (24 before).

## Wider Windows testing: /arch:AVX2, ARM64, Windows bit identity; MSVC /fp:contract fix

Details are in the README section *Supported toolchains*.

* **CI.** The *windows* job gained six builds: on x64 MSVC `/arch:AVX2`, MSVC `/arch:AVX2 /fp:contract` and clang-cl `/arch:AVX2` (Release; the job first checks that the CPU has AVX2), and on ARM64 (`windows-11-arm`, native) MSVC Debug and Release and clang-cl Release. New jobs *bit-identity (windows x64)* and *bit-identity (windows arm64)* (`tests/bitident/windows.py`). CI now has 24 jobs (16 before). The nightly fuzz workflow's MSVC job became a matrix of five Windows builds (double): MSVC (still named *fuzz (msvc)*), MSVC `/arch:AVX2 /fp:contract`, clang-cl `/arch:AVX2`, MSVC ARM64 and clang-cl ARM64; the workflow has 10 jobs (6 before).
* **Results.** All 9 test suites pass in all nine Windows jobs. Fuzz at seed 20261014, 250 samples per function (double): 0 failures in all five Windows builds; every worst error within 1.37× of plain MSVC's (`expint` 9.18 vs 6.68 ε, MSVC `/arch:AVX2 /fp:contract`) and within 1.79× of GCC 14's (`log` 1.62 vs 0.907 ε at z ≈ 1.571 + 0.00057i, in all five Windows builds). This corrects the README's earlier "within 1.5× of GCC 14" for MSVC, which no longer held for `log` after the |z| = 1 change (GCC 14 went to 0.907 there; MSVC stayed at 1.62).
* **Fix: MSVC `/fp:contract`** (no change on GCC, Clang or clang-cl, or on MSVC without `/fp:contract`). MSVC fuses multiply-adds across statements and inlined calls, e.g. the product in `two_prod` with the sum in the `fast_two_sum` that follows, which broke the double-length (`dd`) arithmetic of the special functions. Both headers now compile their own code with `#pragma fp_contract(off)` when `_M_FP_CONTRACT` is defined (MSVC `/fp:contract` and `/fp:fast`) and restore `fp_contract(on)` at their end; template instantiations keep the state of their definition (checked: the suites pass with and without the restore).
  * Before → after, fuzz at seed 20261014, 250 samples per function, double, x64 `/arch:AVX2 /fp:contract` / ARM64 `/fp:contract`: failures 88 / 89 → 0 / 0; `cyl_bessel_j` 1.33e5 / 1.02e5 → 24.7 / 24.7 ε, `sph_bessel` 4.71e4 / 4.71e4 → 8.13 / 8.88, `riemann_zeta` 1.19e3 / 1.19e3 → 4.02 / 4.02, `cyl_bessel_j_int` 1.17e3 / 1.17e3 → 64.1 / 64.1, `cyl_bessel_i` 871 / 968 → 22.8 / 22.8, `tgamma` 60.7 / 60.6 → 2.51 / 2.51, `erfc` 28.8 / 28.8 → 4.39 / 4.39, `airy_bi` 90.4 / 209 → 32.7 / 32.7.
  * Test suites (x64 `/arch:AVX2 /fp:contract`, Release): before, 281 failed checks in 3 suites (`test_special` 24: J_n at 15 for n = 0, 1, 3 and −2; `test_special2` 252: `tgamma`, `erf`, `erfc`, `riemann_zeta`, `cyl_bessel_j` against `std::`, the real-axis mpmath tables, Bessel near zeros; `test_constexpr` 5: the runtime error-free helpers); after, all 9 suites pass.
  * Bit identity: with the pragma, MSVC `/fp:contract` gives the plain MSVC build's bits for the transcendental-free operations (the Linux goldens). Against the pre-constexpr header, which had no pragma, 30 of the 45 hashes differ in that build only; the Windows bit-identity job lists this and does not fail on it.
  * MSVC `/fp:fast` stays unsupported: `test_constexpr`'s `static_assert`s do not compile with it (C2607).
* **Found, no change needed.**
  * MSVC does not contract with `/arch:AVX2` alone: bit-identical to plain MSVC (old and new header, all 45 hashes); with it, `exact_mul` (`log` near |z| = 1) uses `std::fma` (`__AVX2__`).
  * clang-cl contracts by default with `/arch:AVX2` and on ARM64 (30 of 45 hashes differ from the non-contracting builds, as on Linux with `-march=haswell`; identical with `-ffp-contract=off`); old and new header still agree bit for bit, all suites pass, and fuzz is clean.
  * The ARM64 Universal CRT's `hypot` rounds differently from the x64 one: double `abs`/`normalize` hashes `2d7f00354389fa6a`/`64940283057f2787` (new `double@ucrt-arm64` goldens); it differs from a correctly rounded `hypot` in 29,052 of 200,000 double inputs (x64: 32,576) and 370 of 200,000 float inputs (same as x64); a constexpr `abs` differs from the runtime one in 23 of 256 values (x64: 25), 1 ulp each (`test_constexpr`, seed 42). Every other hash matches the glibc goldens on ARM64.
* **Tests.**
  * `QTEST_FP_CONTRACT` (`tests/test_framework.hpp`) decides where the bit-exact runtime checks are skipped: `__FMA__`, or Clang on ARM64/AArch64 (was `__FMA__` only).
  * `test_bitident`: `double@ucrt-arm64` goldens on ARM64; no unreachable code after the FMA skip (MSVC C4702 under `/WX` with `/fp:contract`).
  * `test_constexpr` is built with `-ffp-contract=off` under clang-cl too (`/clang:-ffp-contract=off`), as already under GCC and Clang.
  * The *windows* jobs print `test_constexpr`'s `hypot` statistics after the tests.
* **Not tested on Windows:** sanitizers, benchmarks, MinGW, 32-bit x86, `/arch:AVX512`, ARM64EC.

## log and log10 relatively accurate next to |z| = 1 (changes log, log10 and pow results)

Details are in the README section *log near |z| = 1*.

* **Change.** For a finite q with a nonzero vector part and 0.49 ≤ ‖q‖² ≤ 2.02, ln‖q‖ is now ½·log1p(‖q‖² − 1), with ‖q‖² − 1 summed in double length from the exact squares (`std::fma` where the build has hardware FMA, Dekker's split otherwise and for x87 `long double`) by repeated error-free cascades. Before, it was log(abs(q)), accurate only absolutely there. The same code runs on GCC, Clang, MSVC and clang-cl; elsewhere nothing changed.
* **Results that change** (only for such q; real, zero and non-finite arguments are unchanged):
  * `log`, `log10`: the real part, by at most about ε absolutely; the vector part is bit-identical (200,000 random quaternions with 0.5 ≤ ‖q‖ ≤ 2). 19,980 of the 20,000 `log_unit_check.py` points differ in at least one type.
  * `pow(q, x)`, `pow(q, p)`, `pow(q, complex)`, `pow(complex, p)` with such a base: every component, by about |p|·ε relative (11,373 of 20,000 results differ for p in [−4, 4], 17,986 of 20,000 for |p| from 10 to 1000). `pow(T, p)`, `pow(q, int)` and all other functions are unchanged, and the other rows of the seed-20261014 fuzz report are identical.
  * The bit-identity goldens cover only the transcendental-free operations: unchanged, no hash updated.
* **Accuracy** (worst relative error in ε, before → after):
  * Fuzz, seed 20261014, 250 points per function and type, all four fuzz builds: double `log` 1.47e6 → 0.907, `log10` 6.8e5 → 1.26; long double 1.14e8 → 1.14, 1.78e8 → 1.28. With the new |z|~1 region (real part checked on its own): double 4.5e15 → 0.907 and 4.5e15 → 1.26, long double 9.22e18 → 1.14 and 9.22e18 → 1.36; 137 failures → 0.
  * `fuzz/log_unit_check.py log`, 20,000 points, seed 1, real part: double 3.56e16 → 0.98 (`log`), 3.56e16 → 1.68 (`log10`); long double 9.22e18 → 1.84, 9.22e18 → 2.09; float 8.35e8 → 0.998, 8.35e8 → 1.26. The same with GCC 14 (with and without `-march=haswell`), Clang 19 `-march=haswell` and Clang 18 with libc++.
  * `pow` (`log_unit_check.py pow`, 20,000 points): the normwise worst case is set by the angle p·arg q and stays (p in [−4, 4], seed 5: double 20.3 → 20.3, long double 22.1 → 22.2; |p| from 10 to 1000, seed 6: 5.38e3 → 5.33e3, 4.78e3 → 4.78e3); the error of ‖q^p‖ drops (4.23 → 1.71 and 3.42 → 2.98; 702 → 266 and 713 → 311). Fuzz `pow` stays at 16 (double) and 21.4 ε (long double).
* **Fuzz.** The listed-but-not-failed exemption for `log`/`log10` "complex near 1" is gone (`KNOWN` is empty), so these fail the job again. New region for `log`, `log10` and `pow`: points next to the whole unit circle (|z| = 1 ± 10⁻¹⁷…10⁻², a quarter with x² + y² − 1 cancelling beyond double length), with the real part of `log`/`log10` also checked on its own. Its points come from a generator seeded by the point they replace, so every other function keeps its samples at a given seed. New script `fuzz/log_unit_check.py` with evaluator `fuzz/log_unit_eval.cpp` (not part of CI).
* **Speed** (medians of 5 alternating runs, `-O3`, new/old, GCC 14 / GCC 14 haswell / Clang 19, double): `log` 0.90 / 0.85 / 0.98, `log10` 0.93 / 0.89 / 0.95, `log_unit` 0.97 / 0.85 / 0.88, `log_far` (unchanged path) 0.97 / 1.01 / 1.05, `pow` 0.91 / 0.89 / 0.99, `pow_unit` 1.00 / 0.89 / 0.98; long double 0.85–1.04. The new path skips abs(q) and its three `hypot` calls (GCC 14 double `log` 61.0 → 54.9 ns). Unchanged `exp` read 1.41 with GCC 14 in double, the noise of the shared machine.
* **Tests.** `test_math` checks five points next to ‖q‖ = 1 against mpmath to 4 ε (double, long double); 22 of those checks fail with the old header. New bench cases `log_far`, `log_unit`, `pow_unit`.

## MSVC and clang-cl support (Windows)

Details are in the README section *Supported toolchains*.

* **CI.** New job *windows* on `windows-latest` with CMake's default Visual Studio generator: MSVC 19.51 (toolset v145) in Debug and Release, and clang-cl 22.1 in Release, all with `/std:c++latest`. CI now has 16 jobs (13 before). The nightly fuzz workflow gained a *fuzz (msvc)* job (Release, double).
* **Results.** All 9 test suites (core, math, io, special, special2, annexg, constexpr, bitident, example) pass in all three Windows jobs, including the comparisons with the MSVC STL's `std::` special math functions and every `constexpr` `static_assert`. Fuzz at seed 20261014, 250 samples per function (double): 0 failures; every function within 1.5× of GCC 14's worst error at the same seed (largest ratio `atan`, 1.56 vs 1.07 ε).
* **Source fixes** (no change on GCC or Clang):
  * `two_prod`: the literal `0x1p+16000L` in the x87-only branch is out of range where `long double` is binary64, and MSVC rejects it even in a discarded `if constexpr` branch (C2177, in all 9 targets). It is now `cx::pow2<T>(16000)`.
  * The converting constructor `quaternion<T>(const quaternion<U>&)` decides `explicit` by floating-point conversion rank (`detail::fp_widens`) instead of a requires-expression inside `explicit(...)`. MSVC evaluated that form so that `quaternion<float>` was not constructible from `quaternion<double>` at all (C2607 in `test_core`). Same explicit-ness on every compiler.
  * Vector-part norm on Windows: `detail::vnorm3`, correctly rounded, instead of `hypot(hypot(x, y), z)` from the Universal CRT, whose `hypot` is up to 1 ulp off even when the result is representable. Before → after, MSVC fuzz at seed 20261014: `sinh` at ‖v‖ = 109 63.8 → 1.75 ε (bound 16), `tgamma` at ‖v‖ = 25.6 69.1 → 2.51 ε (bound 64), `riemann_zeta` at ‖v‖ = 237 and 250 482 and 486 → 4.02 ε overall (bound 64), `tgamma_upper` 6.01 → 2.24 ε. The CRT's other functions were not the cause: over 4,000 inputs each (seed 20261014) its worst errors are `sin` 0.66, `cos` 0.74, `sinh` 0.85, `cosh` 0.94, `exp` 0.52, `log` 0.5, `expm1` 2.27 and `log1p` 0.99 ulp (glibc: 0.51, 0.51, 1.26, 0.98, 0.5, 0.5, 0.72, 0.74). `vnorm3` was correctly rounded in all of 3,000 double and 2,979 float test norms (1,000 of them exact oblique-axis norms), at about the speed of the nested glibc `hypot` (1.0–1.35×, noisy). Linux builds keep `std::hypot`.
* **Tests.**
  * `test_bitident`: where `long double` is binary64 its hashes must equal double's; double `abs` and `normalize` (C library `hypot`) have Universal CRT hashes. Before: 16 of 90 checks failed on clang-cl.
  * `test_constexpr`: every difference between a constexpr and a runtime `abs` must be at most 1 ulp (new check); the bound of 16 differences in 256 now applies to glibc only. On the Universal CRT 25 of 256 differ; its `hypot` differs from a correctly rounded one in 370 of 200,000 float and 32,576 of 200,000 double inputs (seed 42).
  * Tests and `fuzz_eval` build with `/W4 /WX` on MSVC. That needed: an `int` cast in `tests/bitident/ops.hpp` (C4805), `if constexpr … else` in two `test_special2` templates (24 × C4702 in Debug), `_CRT_SECURE_NO_WARNINGS` for `fuzz_eval`'s `scanf`, and a `-Wliteral-range` suppression for the x87-only reference 1.01e-3005 under clang-cl.
* **Limits on Windows.** `long double` has double precision. `abs`/`normalize` use the CRT's `hypot` at run time. Sanitizers, benchmarks, `/arch:AVX2` builds, MinGW, 32-bit x86 and ARM64 are not tested. macOS is still on hold.

## FMA audit: near-equal subtractions, elementary-function fuzzing, haswell fuzz builds

After the Hurwitz ζ fix (1adef77), every place that subtracts nearly equal values was checked in FMA and non-FMA builds. Details, the full list of sites and the documented limits are in the README section *FMA audit*.

* **Method.** New script `fuzz/fma_audit.py` runs the same exact inputs through four builds of `fuzz/fuzz_eval.cpp` (GCC 14 and Clang 19 at `-O2`, each without and with `-march=haswell`) and compares each with mpmath. It covers 50 targeted sites (`--sites`) and 64 functions on the fuzz regions. All audit figures: double, worst error in ε, 300 points per site or function, seed 1.
* **Fixed** (before → after):

  | Site | Before | After |
  |---|---|---|
  | `erfc`, 0.4 ≤ \|Re z\| ≤ 1.3 | 19.4 (no FMA), 47.9 (GCC haswell); flagged FMA-sensitive | 8.46, 12.3 (10.6 Clang haswell) |
  | `erf`, tiny vector part | 76.5 | 17.7 |
  | `comp_ellint_1`, k next to ±1 | 7.28e5 | 2.45 |
  | `comp_ellint_3`, real k next to ±1 | 7.6e5 (no FMA), 6.54 / 6.27 (haswell); flagged FMA-sensitive | 6.66 / 6.54 / 6.27 |
  | `ellint_1`, k next to 1, φ next to π/2 | 5.39e11 | 4.68 |
  | `ellint_2`, φ next to π/2 | 8.35 | 7.64 |
  | `ellint_3`, k next to 1, ν sin²φ next to 1 (÷ conditioning in φ) | 1.19e6 | 0.406 |
  | `atanh`, z next to ±1 with \|z\| > 1 | 8.5e4 | 0.471 |

  * `erf`/`erfc`: the double-length Maclaurin series now covers |Re z| < 5/4 with |z| < 4 (was |Re z| < ½), and `erfc` there is 1 − erf from the double-length sum, rounded once. The Laplace continued fraction needed about 700 steps near |z| = ½, and FMA contraction changed how their rounding added up.
  * Elliptic integrals: k′² = (1 − k)(1 + k) instead of 1 − k·k (which is exact only if contracted into `fma(−k, k, 1)`); 1 − k² sin²φ = k′² + (k cos φ)²; the reduced amplitude's sine and cosine are (−1)^m sin φ, (−1)^m cos φ instead of sin, cos of φ − mπ; and in Π(ν, φ, k), 1 − ν sin²φ is (1 − Re ν) + Re ν cos²φ where |cos φ| ≤ |sin φ|, otherwise formed with sin²φ and its product with Re ν carried exactly.
  * `atanh` (and `atan` = −i·atanh(iz)): Kahan's forms ¼·log1p(4x/((1−x)² + y²)) and ½·atan2(2y, (1−x)(1+x) − y²) for |Re z|, |Im z| < 2, instead of (log1p z − log1p(−z))/2 and the detour through 1/z.
  * Over the 64 functions, `erf` (14.6 → 7.66 without FMA, 29.4 → 8.58 with GCC haswell) and `erfc` (22.8 → 7.77, 49.3 → 4.7) were the only ones flagged FMA-sensitive; none is now. Also `comp_ellint_1` 9.46e5 → 3.2, `ellint_1` 5.94e7 → 3.51, `ellint_2` 76 → 56.7, `comp_ellint_2` 38.6 → 36.4, `atanh` 1.05e5 → 1.01.
* **Found by the new elementary-function fuzz regions** (`fuzz/fuzz.py`, 250 points per function and type, seed 20261014), both fixed:
  * Real arguments of `asin`, `acos`, `asinh`, `acosh`, `atanh` and `atan` take the Annex G path, which called the `std::complex` functions through argument-dependent lookup instead of this library's kernels; on libc++ they lost accuracy (50 failures). Clang 18 with libc++, before → after: double `asin` 3.23e3 → 1.87, `asinh` 812 → 1.65, `acosh` 2.29e6 → 1.27, `atanh` 912 → 1.38; long double `asin` 73.7 → 1.68, `asinh` 2.78e3 → 1.56, `acosh` 4.04e7 → 1.4, `atanh` 2.67e3 → 1.42.
  * `atan` with a tiny vector part (every library): double 4.5e15 → 1.07 (1.08 with haswell), long double 9.22e18 → 1.51.
* **Documented limits, not FMA effects:** `log`/`log10` next to |z| = 1 (log(hypot), absolute accuracy: up to 1.47e6 / 6.8e5 ε in double and 1.14e8 / 1.78e8 ε in long double; listed by the fuzz job but not failed); the vector part of `expint`, `cyl_bessel_j` (|z| > 16), `airy_ai` and `polylog` for tiny |v| (absolute, about ε·|f|); `comp_ellint_2` next to k = 1 (46.5 ε before and after, from R_F − (k²/3)R_D).
* **Nightly fuzz.** `fuzz.py` now also samples the 17 elementary functions (16 ε bounds, 64 ε for `pow`), including tiny vector parts and points next to ±1, 0, π/2 and π, and puts a quarter of the points on the oblique axis (2, 3, 6)/7. The workflow gained two FMA builds, GCC 14 and Clang 19 with `-march=haswell` (after a check that the runner has FMA and AVX2), and an informational *fma-audit* job. Locally all four fuzz builds pass with no failures at seed 20261014. A harness bug was fixed on the way: two envelope lines had landed in `reference()`, where they used an undefined variable, so every `expint` point and every Hurwitz ζ point with Re s < 0 raised an exception and was skipped as an mpmath failure. Those points are now checked (the Hurwitz term is |a + n|^−s, with the arguments the right way round).
* **Speed** (medians of 5 alternating runs of the old and new header, `-O3`, new/old, GCC 14 / GCC 14 haswell / Clang 19): `erf` 0.55–0.61, `erf(0.8+1.5i)` 0.28–0.42, `erfc(0.6+0.5i)` 0.08–0.12; `atanh` 0.35–0.45 and `atan` 0.43–0.51; unchanged `erfc` paths 0.97–1.06; the six elliptic cases 0.92–1.11 (shared machine, so within noise).
* **Tests.** All suites pass with GCC 14, Clang 18 and 19 on libstdc++ and libc++ in Debug and Release, with sanitizers, and with GCC 14 and Clang 18 at `-march=haswell`; the haswell bit-identity comparison against 91cad1a is unchanged (45 hashes, `-O2` and `-O3`).
* **New bench cases:** `erf(0.8+1.5i)`, `erfc(0.6+0.5i)`, `erfc(2.5+i)`, `erfc(0.6+5i)`, `comp_ellint_1`, `comp_ellint_2`, `comp_ellint_3(0.6)`, `ellint_1(0.8)`, `ellint_2(0.8)`, `ellint_3(0.6,1.2)`.

## Fix: Hurwitz ζ with FMA contraction (GCC 14 `-march=haswell`)

* **Symptom.** With GCC 14 and `-march=haswell`, `test_special2` failed 13 `hurwitz_zeta` checks:
  * float and double returned `-inf` at s = −1.75 + 120i and −1.25 + 120i (a = 1.25);
  * double errors reached 25–103 |s|ε at |Im s| = 20 and 120.
* **Root cause.** The Hermite-integral path formed (a − it)^{−s} − (a + it)^{−s} as exp(−s·log(a − it)) − exp(−s·log(a + it)) and divided by e^{2πt} − 1 ≈ 2πt. For real `a`, the two logarithms are exact conjugates, and without FMA the products s·log(a ∓ it) round symmetrically, so their difference vanishes like t.
  * GCC contracts the complex product into `fma(si, L, ±sr·θ)`. When si·L lies at or next to a rounding tie, the tiny ±sr·θ (θ ≈ t/a) breaks the tie in opposite directions.
  * The imaginary parts of the two exponents then differ by 1 ulp for every small t (shown with t = 1e-320 … 1e-30: `-0x1.7e9c28928a349p+8` against `…348p+8`).
  * That difference, roughly ε·|a^{−s}|, divided by 2πt overflowed at the double-exponential nodes near t = 1e-320, which gave `-inf`. At larger t it added about ε/t per node, which gave the 25–103 |s|ε errors.
  * Clang contracted differently and happened to pass.
* **Fix (contraction-robust, no flags).** The difference is now e^{−s·c}·(e^{s·d} − e^{−s·d}), with c = log|a + it| and d = i·atan(t/a) for real `a`. For complex `a` (polylog), d = atanh(it/a) = (log1p(w) − log1p(−w))/2, which has no cancellation.
  * It is evaluated as 2·e^X·sinh(Y) when |Re Y| ≤ 1; otherwise the two exponentials are formed separately as before, so nothing can overflow.
  * The small difference is no longer obtained by subtraction, so it no longer depends on how the products are rounded.
* **Errors against mpmath** (161 references in `test_special2`, worst in |s|ε):

  | Build | Before | After |
  |---|---|---|
  | GCC 14 `-march=haswell`, float | inf | 0 |
  | GCC 14 `-march=haswell`, double | inf (2 × `-inf`; others up to 103) | 2.58 |
  | GCC 14 `-march=haswell`, long double | 3.69 | 3.69 |
  | Without FMA, double / long double | 2.58 / 3.69 | 2.58 / 3.69 |

  A random sweep of 150 points (Re s ∈ (−2, 0.5), |Im s| ∈ [8, 150], a ∈ {0.375, …, 4.5}) gives:
  * double: max 1.89 → 1.87 without FMA, and 2.27 → 1.87 with FMA;
  * long double: 2.1 → 1.78.
* **Speed.** The real-`a` path now needs one real log and one atan2 per node instead of two complex logs. Medians of 5 alternating runs, `-O3`:

  | Case | double, new/old | long double, new/old |
  |---|---|---|
  | `hurwitz(−1/2+20i, 2.75)` | 0.88 (GCC), 0.54 (GCC haswell), 0.91 (Clang 19) | 0.91–0.93 |
  | `hurwitz(−7/4+120i, 1.25)` | 0.20–0.23 | 0.14 |
  | `polylog(1/2, −3+i)` (unchanged path) | 0.98–1.03 | 0.88–1.03 |

  The large-|Im s| case was slow before because the cancellation noise made the double-exponential rule refine further.
* **Other checks.**
  * `test_bitident` still passes; the arithmetic is untouched.
  * The haswell CI job no longer skips `special2` for GCC.
  * New bench cases: `hurwitz(...)` ×2 and `polylog(1/2,-3+i)`.

## constexpr without transcendentals (step 1)

* **Now `constexpr` on GCC 14, Clang 18/19, libstdc++ and libc++:**
  * all arithmetic, including the Annex G division (previously `p * inverse(q)` in constant expressions; it now runs the same scaled algorithm as at run time);
  * `abs`, `normalize` and `proj`, joining `norm`, `conj`, `inverse`, `dot`, `cross`, integer `pow` and the comparisons, which already were.
* **New internal helpers** `quat::detail::cx::{isnan, isinf, isfinite, signbit, fabs, copysign, logb, ilogb, scalbn, ldexp, frexp, sqrt, hypot}`.
  * They use `std::bit_cast` for float, double and x87 `long double`, with an exact arithmetic fallback for other formats.
  * They are used only in constant evaluation (`if consteval`); at run time the same `std::` calls as before run.
  * `scalbn`, `logb`, `frexp` and `sqrt` match glibc bit for bit on 200,000 inputs per type.
  * `hypot` is correctly rounded. It differs from glibc's `hypot` in 0.2–0.4% of inputs (1 ulp each), and in all 300 double cases checked with mpmath ours was the closer one.
* **Runtime bit-identity (condition 1).**
  * `test_bitident` compares hashes of 15 operation groups (1,600 input pairs per type, including special values, subnormals and near-overflow values) with golden hashes from commit 91cad1a. It passes on every CI toolchain in Debug, Release and with sanitizers.
  * The new CI job *bit-identity (haswell)* compiles the corpus against the old and the new header with GCC 14 and Clang 18 at `-O2`/`-O3 -march=haswell` and diffs the raw bits; locally all four comparisons are identical.
  * FMA builds have no golden hash, because their contraction pattern changes with the optimization level.
* **No contraction change (condition 2).**
  * No `-ffp-contract=off` or pragma was added to library code. Disabling contraction in the runtime code would itself change results on FMA builds and break condition 1. Only the test target `test_constexpr` is compiled with `-ffp-contract=off`, because it runs the helpers' error-free transformations at run time.
  * Benchmarks: `bench --filter core`, 12 cases × 3 types, median of 5 alternating runs of the old and new header, `-O3`. New/old geometric mean (min–max):

    | Build | Geometric mean | Range |
    |---|---|---|
    | GCC 14 | 0.999 | 0.970–1.047 |
    | GCC 14 `-march=haswell` | 0.990 | 0.911–1.030 |
    | Clang 19 | 0.992 | 0.940–1.048 |
    | Clang 19 `-march=haswell` | 1.010 | 0.957–1.056 |

    The spread is run-to-run noise. Examples, double, ns/call (old → new): `div` 46.9 → 47.4 (GCC) and 43.6 → 44.0 (Clang); `abs` 18.7 → 18.7 (GCC) and 18.9 → 19.0 (Clang); long double `div` 377.7 → 375.6 (GCC).
  * New bench cases: `scalar_div`, `norm`, `normalize`, `proj`, `pow_int`, `pow_neg_int`.
* **Limits.**
  * In FMA builds, a constexpr product or quotient can differ in the last bits from the contracted runtime one.
  * A constexpr `abs` can differ from the runtime one by 1 ulp: the constexpr `hypot` is correctly rounded, glibc's is not always.
  * Operations that would produce a NaN (or, on GCC, overflow) are not constant expressions, so Annex G special-value results are runtime-only.
* **Found while testing (pre-existing, not changed here).** With GCC 14 and `-march=haswell`, `hurwitz_zeta` failed 13 checks in `test_special2`. Fixed since: see "Fix: Hurwitz ζ with FMA contraction" above.
* **Tests.**
  * `test_constexpr`: `static_assert`s for the helpers and the operations.
  * Checks that 256 constexpr quotients per type match the runtime ones on non-FMA builds.
  * The helper sweeps described above.
* **README.** New "constexpr" section with the feasibility study: GCC folding coverage, Clang/libc++ coverage, P0533/P1383 status, special values that never fold, fold-vs-runtime differences of 1,118/16,800 (double) and 1,604/16,800 (long double), 2,047/4,000 contraction differences with `-march=haswell` against 0 with `-ffp-contract=off`, and the 1,500–2,500 line estimate. Also the plan for steps 3 and 4.

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
  * The Airy bound is now 128 ε, because the condition number |z·Ai′(z)/Ai(z)| grows like |z|^{3/2}. The largest error seen was 69.6 ε in long double at z ≈ −18.6 + 5.8i (clang-18/libc++), where mpmath gives a condition number of 86.2 (|z|^{3/2} = 86.3).
  * Over 200 inputs perturbed by up to 10⁻³ around that point, the errors stayed below 1.03 × condition × ε, with a median of 0.3. This held for double and long double, and for libstdc++ and libc++. The same point gives 8.4 ε with gcc-14/libstdc++. The size of the error follows the conditioning rather than any one toolchain.
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
