# Changelog

Accuracy figures are worst-case errors against mpmath on random exact inputs (300 points per region unless noted), double / x87 long double.

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
