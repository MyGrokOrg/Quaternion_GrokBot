# Changelog

Accuracy figures are worst-case errors against mpmath on random exact inputs (300 points per region unless noted), double / x87 long double.

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
