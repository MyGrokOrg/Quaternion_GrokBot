# Quaternion_GrokBot

Header-only **C++23 quaternions** (`a + bi + cj + dk`) whose API follows `std::complex`.

```cpp
#include <quaternion.hpp>          // everything; or <quaternion/quaternion.hpp> for the core only
using namespace quat::literals;

quat::quaternion<double> p = 1.0 + 2.0_i + 3.0_j + 4.0_k;
auto q = quat::quaternion(0.5, -1.0, 0.25, 2.0);   // CTAD -> quaternion<double>

auto r = p * q;                   // Hamilton product (p*q != q*p in general)
std::println("{:.3f}", exp(p));   // (1.694,-0.790,-1.184,-1.579)
std::cout << sqrt(p) << '\n';     // (1.79961,0.555675,0.833512,1.11135)
auto [w, x, y, z] = p;            // structured bindings
```

## Requirements and build

* A C++23 compiler and standard library with `<format>`. It is tested with GCC 14, and with Clang 18/19 on both libstdc++ and libc++. The comparisons against the `std::` special math functions only run where the library provides them (libstdc++).
* CMake 3.20 or newer (optional; the library is a single header).

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
cmake --install build --prefix /your/prefix     # optional
```

Consume it with CMake:

```cmake
find_package(Quaternion 1.0 REQUIRED)            # after install
# or: add_subdirectory(Quaternion_GrokBot)
target_link_libraries(app PRIVATE Quaternion::quaternion)
```

## API reference

Everything lives in namespace `quat`. The table follows the layout of `<complex>`.

### `template<std::floating_point T> class quaternion`

Like `std::complex`, `T` is `float`, `double` or `long double`. The layout is four `T`s (trivially copyable).

| Member | Description |
|---|---|
| `value_type` | `T` |
| `quaternion(T a = 0, T b = 0, T c = 0, T d = 0)` | `a + bi + cj + dk` |
| `quaternion(std::complex<T> z1, std::complex<T> z2 = {})` | Cayley–Dickson pair `z1 + z2·j` = `(z1.re, z1.im, z2.re, z2.im)`. Implicit, so a complex number converts to a quaternion. |
| `quaternion(const quaternion<U>&)` | Converts between precisions. It is `explicit` when the conversion narrows, as `std::complex` does. |
| `real()`, `real(T)` | Read or set the real (scalar) part. |
| `unreal()` | Vector part as a pure quaternion `(0, b, c, d)`. |
| `R_component_1()`…`R_component_4()` | The components `a`, `b`, `c`, `d` (Boost-compatible names). |
| `C_component_1()`, `C_component_2()` | `complex(a, b)` and `complex(c, d)`. |
| `w() x() y() z()`, and setters `w(T)` … `z(T)` | Short names for the components. |
| `get<I>()` | Reference to component `I` (0–3). |
| `= += -= *= /=` | With `T`, `std::complex<T>`, or `quaternion<U>`. `q *= r` means `q = q*r` and `q /= r` means `q = q*r⁻¹` (right multiplication). |
| unary `+`, `-` | |

All members are `constexpr` and `noexcept`. Accessors are `[[nodiscard]]`.

### Non-member operators

* `+ - * /` for every combination of quaternion, `T` and `std::complex<T>` (they are hidden friends, so `q * 2` works).
* `p / q` is `p·conj(q)/‖q‖²` with `q` first scaled by a power of two (as C's `_Cdivd`), so it does not overflow or underflow for large or tiny `q`. In a constant expression it is `p * inverse(q)`. `s / q` is `s * inverse(q)`.
* `*` and `/` recover infinities as C11 Annex G (G.5.1) does; see [Infinities, NaNs and signed zeros](#infinities-nans-and-signed-zeros-c-annex-g).
* `==`, and `!=` through C++20 rewriting, against quaternion, `T` and `std::complex<T>` in either order.
* `operator<<` prints `(a,b,c,d)`. It honors the stream's flags, precision and locale, and the width applies to the whole text, as with `std::complex`.
* `operator>>` reads `a`, `(a)`, `(a,b)`, `(a,b,c)` or `(a,b,c,d)`. On malformed input it sets `failbit` and leaves the target unchanged.

### Value functions

| Function | Result |
|---|---|
| `real(q)`, `unreal(q)` | Scalar part, and the vector part as a pure quaternion. |
| `abs(q)` | `‖q‖`, computed with `hypot` so it doesn't overflow. |
| `arg(q)` | `atan2(‖v‖, a)` ∈ [0, π]: the angle from the positive real axis. |
| `norm(q)` | `a²+b²+c²+d²` (`constexpr`; squared magnitude, as in `std::norm`). |
| `conj(q)` | `(a, −b, −c, −d)` (`constexpr`). |
| `proj(q)` | `q`, or `(+∞, ±0, ±0, ±0)` if any component is infinite (as `std::proj`). |
| `polar(rho, theta = 0)` | `rho(cos θ + i sin θ)`, identical to `std::polar` embedded on the i axis. |
| `polar(rho, theta, axis)` | `rho(cos θ + u sin θ)` with `u` = the unit vector part of `axis`. This inverts `abs`/`arg`: `polar(abs(q), arg(q), q) == q`. |
| `inverse(q)` | `conj(q)/norm(q)` (`constexpr`). |
| `dot(p, q)` | 4-D Euclidean inner product (`constexpr`). |
| `cross(p, q)` | Cross product of the vector parts, as a pure quaternion (`constexpr`). |
| `normalize(q)` | `q/‖q‖`, or `q` when `q` is zero. |

### Transcendental functions

Write `q = a + u·θ`, where `u` is the unit vector of `q` and θ = ‖v‖. Every function below keeps `q`'s axis `u`, and on that axis it agrees exactly with the `std::complex` function. Real arguments use the `i` axis, so `log(-1) = πi` and `sqrt(-4) = 2i`. The sign of a zero `i` component selects the side of a branch cut exactly as for `std::complex`: `log(quaternion(-1.0, -0.0, 0, 0)) = −πi`.

| Function | Definition |
|---|---|
| `exp(q)` | `eᵃ (cos θ + u sin θ)` |
| `log(q)` | `ln‖q‖ + u·arg(q)` (principal branch) |
| `log10(q)` | `log(q)/ln 10` |
| `sqrt(q)` | Principal root (real part ≥ 0) |
| `pow(q, int n)` | Repeated squaring (`constexpr`). Negative `n` gives the inverse. |
| `pow(q, T x)` | `exp(x·log q)` |
| `pow(q, p)` | `exp(p·log q)`, the same operand order as `std::pow` for complex. Note that this is not `exp(log q·p)`. |
| `pow(T x, p)`, `pow(q, complex)`, `pow(complex, p)` | The same definition after converting the arguments. |
| `sin cos tan` | `sin(q) = sin a cosh θ + u cos a sinh θ`, and so on. `tan(q) = −u·tanh(u q)`, with Kahan's overflow-free complex tanh. |
| `sinh cosh tanh` | `sinh(q) = sinh a cos θ + u cosh a sin θ`, and so on. `tanh` uses Kahan's algorithm (no overflow for large `a` or θ). |

### Infinities, NaNs and signed zeros (C Annex G)

The elementary functions (`exp log log10 sqrt pow sin cos tan sinh cosh tanh asin acos atan asinh acosh atanh`), the arithmetic operators, `abs`, `arg`, `proj` and `polar` follow C11 Annex G (G.5–G.6) and C23. On the `i` axis they reproduce the Annex G special-value tables of `std::complex` (`cexp`, `clog`, `csqrt`, `csinh`, `ccosh`, `ctanh`, `casinh`, `cacosh`, `catanh`, `cacos`, and through the Annex G identities `csin`, `ccos`, `ctan`, `casin`, `catan`). The library uses its own kernels for these, not the C library's, so the results are the same on libstdc++ and libc++.

A quaternion `q = a + v` is mapped to the complex number `z = a + i·|v|` and back. The vector part plays the imaginary part:

* **Vector part zero** (all of `b, c, d` are ±0). Then `z = a + i·b`, where `b` is the signed `i` component. So `(a, ±0, ±0, ±0)` is the `std::complex` argument `(a, ±0)`, and the result `R + iI` is returned as `(R, I, c′, d′)`. The `j` and `k` zeros keep their signs, and both flip when the sign of `I` differs from the sign of `b`: they behave like an infinitesimal vector parallel to `i`. For example, `sqrt(-4, -0, +0, -0) = (+0, −2, +0, −0)` (no flip), and `acos(2, +0, +0, −0) = (+0, −1.3169…, −0, +0)` (flipped, since Im is negative while b = +0).
* **Nonzero vector part.** Then `z = a + i·|v|` (with `|v| ≥ 0`) and the result is `R + u·I`, where `u = v/|v|`. If `|v|` is infinite, `u` has components `±1/√m` on its `m` infinite components and signed zeros elsewhere. If `v` has a NaN and no infinity, the NaN components receive `I` and the others become signed zeros.
* A zero component of `u` always yields a signed zero, whose sign is the product of the signs of `I` and the component (just the component's sign if `I` is NaN), and never `0·∞ = NaN`. So `exp(−∞, 0, 3, 0) = (−0, +0, +0, +0)` (cexp(−∞ + 3i) = +0·cis 3), and at the pole `atan(0, 0, 1, 0) = (0, +0, +∞, +0)`.
* Annex G entries whose sign is "unspecified" carry that freedom over to the components of `u`.

**Multiplication and division** use the recovery of C11 G.5.1 (`_Cmultd`, `_Cdivd`), generalized to four components. It applies only when the naive result is NaN in all four components:

* `∞·x` with `x` nonzero, or a NaN in the other operand, gives an infinity. The infinite operand is "boxed" (infinite components to ±1, the others to ±0), the product is recomputed and then scaled by ∞.
* `x/0` (with `x` not all NaN) is `∞·x`.
* `∞/finite` is infinite, and `finite/∞` is a (signed) zero.

**Other functions.** `abs` is a nested `hypot`, so it is `+∞` if any component is infinite, even with NaNs. `arg` is `atan2(|v|, a)` ∈ [0, π]; it is never negative, because `|v| ≥ 0`. `proj` is as `std::proj`. `polar(∞, 0)` is `(∞, +0, 0, 0)`, not `∞·0 = NaN`. `pow(q, p) = exp(p·log q)`, and `pow(0, p)` is `1` for `p == 0` and `0` otherwise.

`tests/test_annexg.cpp` checks every Annex G table entry (G.6.1.1–G.6.3.2, with the values ±0, ±½, ±1, ±1.5, ±2, ±∞ and NaN) for `float`, `double` and `long double`. It does this on the `i` axis with both signs of the `j`/`k` zeros, and lifted to the axes `i`, `j`, `−k`, `(1,−1,0)/√2` and `(1,1,−1)/√3` and their negatives. Every sign is checked, except where Annex G leaves it unspecified. `ctanh(±0 ± i∞)` may be `NaN + iNaN` (C11) or `±0 + iNaN` (C23); the library returns the C23 value.

**`-ffast-math` caveat.** `-ffast-math` (and `-ffinite-math-only`, `-fno-signed-zeros`, `/fp:fast`) let the compiler assume that no infinities, NaNs or signed zeros occur. It may then delete the `isnan`/`isinf` tests and the sign handling above. Annex G behavior is only guaranteed without these flags. Finite results are unaffected apart from the usual fast-math reassociation.

### Special functions (`<quaternion/special_functions.hpp>`, included by `<quaternion.hpp>`)

Each function extends an analytic function `f` by **axis reduction**. Write `q = a + |v|·u`, where `u` is the unit vector part. Then

    f(q) = Re f(a + i|v|) + u · Im f(a + i|v|)

So `f(q)` lies in the same complex plane as `q` and commutes with it. On the `i` axis it agrees exactly with the complex function. When `f` is real on the real axis, `f(conj q) = conj f(q)`, and a real argument gives a real result.

| Function | Notes |
|---|---|
| `asin acos atan asinh acosh atanh` | Branch cuts are those of `std::complex` (C99 Annex G). `asin` and `acos` use the algorithm of Hull, Fairgrieve and Tang (1997), with `asinh z = −i·asin(iz)` and `acosh` from `acos`. `atanh` and `atan` have their own kernels: `(log1p z − log1p(−z))/2`, with `atanh(1/z) ± iπ/2` for `|z| > 1`, and `atan z = −i·atanh(iz)`. They are accurate to a few ulp on every standard library. |
| `expm1(q)`, `log1p(q)` | Keep full relative accuracy for small `q`. |
| `cbrt(q)` | Principal cube root (`cbrt(-8) = 1 + √3 i`), like `std::pow(z, 1/3)`. For real `x ≥ 0` it matches `std::cbrt`. |
| `tgamma(q)` | Γ = exp(log Γ), with log Γ carried in double-length (hi + lo) arithmetic, so the size of log Γ costs no accuracy. Uses reflection `Γ(z)Γ(1−z) = π/sin πz` for `Re < ½`, with sin πz kept in split exponential form so it cannot overflow. Gives ±∞ at 0 and NaN at negative integers, like `<cmath>`. |
| `lgamma(q)` | log Γ on the analytic-continuation branch (mpmath's `loggamma`). The real part is `log|Γ|`, which equals `std::lgamma` on the real axis. For `Re z < −4000` it uses reflection, then restores the continuous branch by choosing the multiple of 2πi from Stirling's formula. |
| `digamma(q)` | ψ = Γ′/Γ. Uses the asymptotic series after recurrence, and reflection for `Re < ½`. |
| `beta(q, b)`, `beta(b, q)` | B(q, b) = Γ(q)Γ(b)/Γ(q+b) for real `b`. Two quaternions on different axes have no well-defined beta. |
| `erf(q)`, `erfc(q)` | Maclaurin series for `|Re| < ½`, summed in double-length arithmetic. `exp(−z²)` is formed with `z²` split exactly by FMA. Otherwise the Laplace continued fraction for `erfc` (modified Lentz), with `erfc(−z) = 2 − erfc(z)`. |
| `riemann_zeta(q)` | ζ. On `Re ≥ ½`, uses Borwein's η algorithm for `|Im s| ≤ 1` and Euler–Maclaurin summation otherwise. The Euler–Maclaurin sum takes N ≈ 0.6·|s| + 20 terms plus up to 20 Bernoulli corrections, so there is no cap on `|Im s|`; the cost is O(|s|). On `Re < ½` it uses the functional equation. The phases `Im(s)·log k`, the partial sum and the functional-equation exponent are carried in double-length arithmetic. ζ(1) = ∞, ζ(0) = −½, and the trivial zeros are exact. |
| `cyl_bessel_j(n, q)` | J_n for integer `n` (an entire function). Uses the power series for `|z| ≤ 1`. Otherwise it uses Miller's backward recurrence, normalized by `1 = J₀ + 2ΣJ₂ₖ` or `cos z = J₀ + 2Σ(−1)ᵏJ₂ₖ`, whichever cancels less. |

**Branch cuts.** For a nonzero vector part, the complex kernel is only evaluated on the closed upper half-plane, since `|v| ≥ 0`. On a branch cut along the real axis, the result is the limit from above (imaginary part +0), which is what `std::complex` gives for `(x, +0)`. For `asin acos atan asinh acosh atanh`, a real argument with a −0 `i` component gives the lower side, as `std::complex` does for `(x, −0)`; the other special functions always use the upper side. For example, `asin(2) = π/2 + 1.3169…i`, `lgamma(−2.5) = −0.0562… − 3πi`, and `log1p(−3) = log 2 + πi`. Any imaginary part of a real argument goes on the `i` axis, as with `sqrt` and `log`.

**Accuracy.** The kernels run in at least double precision, so `float` results are close to correctly rounded. Measured against mpmath at 2,800 random points (|Re| ≤ 40, |Im| ≤ 60), the worst relative error was:

| | double | long double (x87) |
|---|---|---|
| tgamma, lgamma, digamma | 2e-14 | 1.5e-16 (close to poles) |
| erf, erfc | 6e-15 | 5e-18 |
| riemann_zeta | 6e-14 | 3e-17 |
| cyl_bessel_j | 1.2e-13 (near real zeros) | 2.5e-17 |

Typical errors are a few ulp. In the large-argument regions, again for exact inputs, these are the worst errors measured against mpmath (300 random points each):

| | double | long double (x87) | region |
|---|---|---|---|
| tgamma | 3.4 ulp | 3.9 ulp | \|Re z\| ≤ 160, 0 ≤ Im z ≤ 300 |
| lgamma | 2.1 ulp | 2.5 ulp | −10⁶ ≤ Re z ≤ −4500, Im z ≤ 10⁵ |
| erf, erfc | 13 ulp | 18 ulp | \|Re z\|, \|Im z\| ≤ 26 |
| asin, acos, asinh, acosh | 2.5 ulp | 3.2 ulp (both libstdc++ and libc++) | 10⁻²⁰⁰ to 10²⁰⁰, next to ±1 |
| riemann_zeta | 19 ulp (≤ 4 ulp for Im s ≥ 1500) | 14 ulp (≤ 4 ulp for Im s ≥ 1500) | −2 ≤ Re s ≤ 3, 1 ≤ Im s ≤ 120,000 |

Some limits are inherent:

* **Γ and log Γ conditioning.** Rounding the input z by ε changes Γ(z) by about `|z·ψ(z)|·ε ≈ |z log z|·ε`. The library adds only a few ulp to that, but no algorithm can remove it for inexact inputs.
* **erf conditioning.** The input sensitivity is about `2|z|²·ε`.
* **ζ conditioning.** The input sensitivity is `|s·ζ′(s)/ζ(s)|·ε`, which is about 7000·ε at |Im s| ≈ 1300 on the critical strip and grows with |Im s|. For exact inputs the library stays within a few ulp. The cost grows linearly in |Im s|: about 1 s at |Im s| = 10⁶ in double, and about 2 s at 10⁵ in x87 `long double`, because `fmal` is done in software.
* J_n for `|z|` beyond about 10⁵: the cost grows like `O(|z|)`.

#### More special functions

These also use axis reduction. Real parameters such as orders, degrees, moduli and `a`, `b` stay real; the argument shown as `q` is the quaternion. Integer orders and degrees are `unsigned`. Floating orders accept any floating type.

| Function | Quaternion argument | Notes |
|---|---|---|
| `cyl_bessel_j(ν, q)`, `cyl_neumann(ν, q)` | `q` | J_ν, Y_ν for real ν (negative allowed), principal branches with the cut on (−∞, 0]. `cyl_bessel_j(int n, q)` still uses the integer-order kernel. |
| `cyl_bessel_i(ν, q)`, `cyl_bessel_k(ν, q)` | `q` | I_ν, K_ν. |
| `cyl_hankel_1(ν, q)`, `cyl_hankel_2(ν, q)` | `q` | H⁽¹⁾ = J + iY and H⁽²⁾ = J − iY, where *i* is the axis of `q`; they are never real. |
| `sph_bessel(n, q)`, `sph_neumann(n, q)`, `sph_hankel_1/2(n, q)` | `q` | j_n, y_n, h_n = √(π/2q) × the order n+½ function. |
| `airy_ai(q)`, `airy_bi(q)` | `q` | Entire. |
| `legendre(l, q)`, `assoc_legendre(l, m, q)` | `q` | P_l and P_l^m = (1−q²)^{m/2} d^m P_l/dq^m with no Condon–Shortley phase, which is the `std::assoc_legendre` convention. Odd m uses the principal square root. |
| `sph_legendre(l, m, θ)` | `θ` | Y_l^m(θ, 0) including the (−1)^m phase, as `std::sph_legendre`. |
| `hermite(n, q)`, `laguerre(n, q)`, `assoc_laguerre(n, m, q)` | `q` | Physicists' H_n and generalized L_n^m. |
| `chebyshev_t(n, q)`, `chebyshev_u(n, q)` | `q` | T_n, U_n. |
| `comp_ellint_1(k)`, `comp_ellint_2(k)` | modulus `k` | K(k), E(k), with the modulus convention of `<cmath>` (m = k²). |
| `comp_ellint_3(k, ν)` | characteristic `ν` | Π(ν, k) = ∫₀^{π/2} dθ / ((1 − ν sin²θ)√(1 − k² sin²θ)), the `<cmath>` convention. Real `k`. |
| `ellint_1(k, φ)`, `ellint_2(k, φ)` | amplitude `φ` | F(φ, k), E(φ, k). Real `k`. The amplitude is reduced modulo π. |
| `ellint_3(k, ν, φ)` | characteristic `ν` | Π(ν, φ, k). Real `k` and `φ`. |
| `expint(q)` | `q` | Ei. On the real axis it is the real principal value, as `std::expint`; elsewhere Ei(z) = −E₁(−z) ± iπ. |
| `expint(n, q)` | `q` | E_n (E₁ is the exponential integral E1). |
| `sinint(q)`, `cosint(q)` | `q` | Si (entire), Ci (cut on (−∞, 0]). |
| `fresnel_s(q)`, `fresnel_c(q)` | `q` | S, C with the kernel sin(πt²/2), cos(πt²/2). |
| `polygamma(m, q)`, `trigamma(q)` | `q` | ψ^(m). |
| `hurwitz_zeta(s, a)` | `s` | ζ(s, a) for real `a`. |
| `polylog(s, q)` | `q` | Li_s for real `s` (integer or not), with the cut on [1, ∞). |
| `lambert_w(q)`, `lambert_w(k, q)` | `q` | W_k for any integer branch k, using the branch conventions of mpmath and Corless et al. W₀ is real on [−1/e, ∞) and W₋₁ is real on [−1/e, 0). |
| `tgamma_upper(a, q)`, `tgamma_lower(a, q)` | `q` | Γ(a, q) for any real a (including a ≤ 0), and γ(a, q) for a not a non-positive integer. |
| `gamma_p(a, q)`, `gamma_q(a, q)` | `q` | Regularized P and Q for a > 0. |
| `ibeta(a, b, q)`, `ibetac(a, b, q)` | `q` | Regularized I_q(a, b) and 1 − I_q(a, b) for a, b > 0. |

**Methods.**

* **K_ν** comes from its Laplace-type integral (DLMF 10.32.8), rotated so the integrand does not oscillate. The integral is evaluated by double-exponential quadrature.
* **H⁽¹⁾_ν** is K_ν(−iz) up to a factor.
* **J_ν** comes from the Wronskian together with the continued fraction for J_{ν+1}/J_ν. I_ν is computed the same way with I_{ν+1}/I_ν.
* **Other quadrants and negative orders** use the analytic-continuation formulas.
* **Airy functions** are computed from K_{1/3}/I_{1/3}, or from J_{1/3}/Y_{1/3} when |ph z| > 2π/3.
* **Elliptic integrals** use Carlson's R_F, R_D and R_J with complex arguments.
* **E₁** uses the power series or the continued fraction. Si, Ci and Fresnel use the power series for small |z|, and otherwise E₁ or erf.
* **ψ^(m)** uses recurrence and then the asymptotic series.
* **Hurwitz ζ** picks one of three methods:
  * Hurwitz's Fourier formula (DLMF 25.11.9) for real `a`, Re s ≤ −2 and |Im s| > 8.
  * Euler–Maclaurin for Re s ≥ ½ and |Im s| ≤ 15.
  * Hermite's integral (DLMF 25.11.29) everywhere else. Before integrating, `a` is shifted by the recurrence so that Re a ≈ max(1, 1.3·|Im s|/2π). That makes the factor e^{|Im s|·atan(t/a)} small next to e^{−2πt}, so the integral does not cancel. The integrand is formed in logarithmic form so it cannot overflow for large |Im s|.
* **Li_s** uses the power series, or the expansion in log z, or inversion. The inversion uses Bernoulli polynomials for integer s and Jonquière's formula otherwise.
* **Lambert W** uses Halley iteration. The start is the branch-point series or the asymptotic expansion.
  * Near −1/e (|W + 1| < ½), the result is refined by Newton iteration on d = W + 1 against the residual e⁻¹g(d) − (z + 1/e), where g(d) = (d − 1)eᵈ + 1 is summed as a series. This form has no cancellation.
  * z + 1/e is formed with 1/e split into two floating-point parts, so it is exact up to one rounding.
* **Incomplete gamma** uses the series or the Legendre continued fraction, and recurrence for a ≤ 0. Incomplete beta uses its continued fraction.
* **Polynomials** are evaluated by their three-term recurrences in quaternion arithmetic, so they are exact polynomials in `q`.

**Accuracy** was measured against mpmath at about 300 random points per function. The table gives the worst relative error.

| | double | long double (x87) | sampled region |
|---|---|---|---|
| J, Y, I, K, H⁽¹⁾, H⁽²⁾ (−12 ≤ ν ≤ 12) | 9e-14 (J near a real zero), otherwise ≤ 2e-14 | 1.2e-17 | \|Re z\|, Im z ≤ 30 |
| Ai | 3e-14 | 3e-17 | ≤ 15 |
| Bi | 5e-12 (near zeros on the negative axis) | 7e-17 | ≤ 15 |
| K, E, F, E(φ), Π | 3e-15 | 1e-18 | ≤ 3 to 6 |
| E₁, Ei, Si, Ci | 6e-15 | 5e-18 | ≤ 30 to 40 |
| Fresnel S, C | 3e-14 | 1e-17 | ≤ 5 |
| ψ^(m), m ≤ 6 | 1.4e-12 (next to poles) | 8e-18 | ≤ 20 |
| Hurwitz ζ, −25 ≤ Re s ≤ 12 | about 50·\|s\|·ε, typically ≤ 10·\|s\|·ε (2e-13 next to a zero in `a`, within that point's conditioning) | same in units of ε | \|Im s\| ≤ 400, 0.05 ≤ a ≤ 5 |
| Li_s, −3 ≤ s ≤ 6 | 1e-12 | 2e-16 | \|z\| ≤ 6 |
| Lambert W, branches −3 … 2 | 2e-16 | 1e-19 | \|z\| ≤ 50 |
| Lambert W next to −1/e (\|z + 1/e\| from 1e-19 to 0.3; branches 0, ±1) | 3 ulp | 3 ulp | |
| Li_s, non-integer 3 < s ≤ 25, e^5.5 < \|z\| < e^40 | 3e-14 | 1.2e-17 | |
| γ(a, z), −4 ≤ a ≤ 8 | 8e-13 | 7e-17 | ≤ 20 |
| Γ(a, z), −8 ≤ a ≤ 20, all arg z (including the negative real axis and a at or near 0, −1, …) | 18 ulp (up to 43 where the input conditioning \|z^a e^{−z}/Γ(a, z)\| is 16–200) | 21 ulp | \|z\| ≤ 100 |
| I_x(a, b), 0.2 ≤ a, b ≤ 6 | 6e-15 | 2e-18 | \|x\| ≤ 1.5 |

Accuracy degrades in these cases:

* **Hurwitz ζ conditioning.** The condition number of ζ(s, a) with respect to s grows like |s|·log, so rounding s alone limits the relative accuracy to about |s|·ε. The library stays within a small multiple of that, up to |Im s| = 400 (the largest tested). Near a zero of ζ(s, ·), the condition number with respect to `a` is |a·s·ζ(s+1, a)/ζ(s, a)|, which can reach 10³–10⁴. That is why the error there can reach about 2e-13; these cases are input conditioning, not algorithm error.
* **Complex `a` (internal).** Li_s calls ζ(1−s, a) with complex `a`, which is accurate there because Im s = 0. With complex `a`, Re s < −6 and |Im s| > 8, Hermite's integral cancels by up to about 10³·|s|·ε. With complex `a`, Re s > 3 and small |Im s|, neither method is reliable. No public function reaches either case.
* **Bessel functions.** The quadrature behind K_ν is adaptive, and the J/I continued fractions need O(|z|) terms, so large |z| or |ν| is slow but still accurate.
* **Lambert W at −1/e.** W is computed to a few ulp of the exact value for the given floating-point input. Because W ≈ −1 + √(2e(z + 1/e)), an input rounding error δ moves the true W by about √(2e)·δ/(2√|z + 1/e|). That error comes from the input, not from this library. Also, a quaternion argument q = a + |v|u can only reach one side of the cut along (−∞, −1/e], namely the |v| ≥ 0 side. So branches k and −k are related by conjugation when the axis is reversed.
* **Γ(a, z).** The Legendre continued fraction is evaluated backwards, from a term count that the modified Lentz algorithm finds first. The forward Lentz product accumulated tens of ulp. The factor z^a e^{−z} is formed with a double-length exponent. The fraction is used unless z is within about 3° of the negative real axis (1 + cos arg z ≤ 0.004), or |z| ≤ 2 (|z| ≤ 0.8 for arg z ≤ 0.6π). Elsewhere a power series is used:
  * For a ≤ ½, the series is taken around the nearest pole a = −n + δ. The pole of Γ(a) and the k = n term are combined analytically as `(−1)ⁿ/n!·[expm1(log Γ(1+δ) − Σ log1p(−δ/m)) − expm1(δ log z)]/δ`, so nothing cancels at or near the non-positive integers.
  * For a > ½, both the Kummer series and the alternating series are evaluated, and the one with the smaller estimated rounding error is kept. If both exceed 16 ε and the fraction is reliable (|z| ≥ a/2), the fraction is used instead.
  * The remaining errors stay within the input conditioning |z^a e^{−z}/Γ(a, z)|.

**Not provided.** Elliptic integrals with two quaternion parameters and Hurwitz ζ with a quaternion `a` are not provided, since axis reduction needs a single quaternion argument.

### Tuple protocol, formatting and literals

* `std::tuple_size<quaternion<T>>` is 4, and `std::tuple_element<I, quaternion<T>>` is `T`. `quat::get<I>(q)` is found by ADL, so structured bindings work, and so does `using std::get; get<0>(q);`.
* `std::formatter<quaternion<T>, CharT>` takes the format spec of `T` and applies it to each component, so `std::format("{:.2f}", q)` gives `(1.00,2.00,3.00,4.00)`. Wide strings work too.
* The inline namespace `quat::literals` (and `quat::literals::quaternion_literals`) defines:
  * `_i`, `_j`, `_k` for `quaternion<double>`
  * `_if`, `_jf`, `_kf` for `float`
  * `_il`, `_jl`, `_kl` for `long double`

  Each accepts both integer and floating literals, as in `2_i` and `3.5_k`.

## Rotations

A unit quaternion `r = polar(1, θ/2, axis)` rotates a vector `v` (stored as a pure quaternion) by θ around `axis`:

```cpp
using Q = quat::quaternion<double>;
Q r = quat::polar(1.0, std::numbers::pi / 4, 1.0_k);   // 90° about z
Q v = r * 1.0_i * conj(r);                              // ≈ (0,0,1,0)
```

## Notes and limitations

* Mixed-precision binary operators such as `quaternion<float> * quaternion<double>` are not provided, which matches `std::complex`. Convert one operand explicitly. The compound assignments do accept other precisions.
* The standard forbids adding overloads to namespace `std`, so `std::get<I>(q)` written fully qualified does not compile. Use `quat::get`, ADL, or structured bindings.
* The special functions lift a complex kernel in one argument. Every other parameter (an order, a degree, the modulus `k`, the Hurwitz `a`, the incomplete gamma/beta `a` and `b`) must be real.
* **libc++ and `long double`.** libc++ 18/19 gives only double-precision `std::numbers` constants for `long double`, and it divides `complex<long double>` values in double precision. Its `long double` complex `asin`, `acos`, `asinh` and `acosh` are good to only about 1e-16, and they overflow for huge arguments. This library therefore carries its own constants, its own complex division (Smith's algorithm), and its own inverse trigonometric and hyperbolic kernels. All of them are a few ulp on every library.
* The transcendental functions are not `constexpr`, because `<cmath>` is not `constexpr` until C++26. All arithmetic, `norm`, `conj`, `inverse`, `dot`, `cross` and integer `pow` are `constexpr`.
* Annex G special values cover the elementary functions and arithmetic, but not the other special functions (gamma, zeta, Bessel, …). With non-finite input these return whatever their kernels produce, usually NaN.
* For large ‖v‖, `exp`, `sin`, `cos`, `sinh`, `cosh` and `tanh` inherit the conditioning of the rounded θ = ‖v‖: a relative error of ε·θ in `cos θ`/`sin θ`. These errors are unchanged by the Annex G work.
* In a `std::format` spec, the width applies to each component, not to the whole quaternion.

## Example

See [`examples/basic.cpp`](examples/basic.cpp), which is built as `example_basic` and run as a test.

## License

MIT. See [LICENSE](LICENSE).
