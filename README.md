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

* A C++23 compiler and standard library with `<format>` (tested with GCC 14, and Clang 18/19 with libstdc++ or libc++).
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
* `p / q` is defined as `p * inverse(q)`. `s / q` is `s * inverse(q)`.
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

Write `q = a + u·θ`, where `u` is the unit vector of `q` and θ = ‖v‖. Every function below keeps `q`'s axis `u`, and on that axis it agrees exactly with the `std::complex` function. Real arguments use the `i` axis, so `log(-1) = πi` and `sqrt(-4) = 2i`.

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
| `sin cos tan` | `sin(q) = sin a cosh θ + u cos a sinh θ`, and so on. `tan = sin·cos⁻¹`. |
| `sinh cosh tanh` | `sinh(q) = sinh a cos θ + u cosh a sin θ`, and so on. |

### Special functions (`<quaternion/special_functions.hpp>`, included by `<quaternion.hpp>`)

Each function extends an analytic function `f` by **axis reduction**. Write `q = a + |v|·u`, where `u` is the unit vector part. Then

    f(q) = Re f(a + i|v|) + u · Im f(a + i|v|)

So `f(q)` lies in the same complex plane as `q` and commutes with it. On the `i` axis it agrees exactly with the complex function. When `f` is real on the real axis, `f(conj q) = conj f(q)`, and a real argument gives a real result.

| Function | Notes |
|---|---|
| `asin acos atan asinh acosh atanh` | Use the `std::complex` kernels and branch cuts. |
| `expm1(q)`, `log1p(q)` | Keep full relative accuracy for small `q`. |
| `cbrt(q)` | Principal cube root (`cbrt(-8) = 1 + √3 i`), like `std::pow(z, 1/3)`. For real `x ≥ 0` it matches `std::cbrt`. |
| `tgamma(q)` | Γ. Uses reflection `Γ(z)Γ(1−z) = π/sin πz` for `Re < ½`. Gives ±∞ at 0 and NaN at negative integers, like `<cmath>`. |
| `lgamma(q)` | log Γ on the analytic-continuation branch (mpmath's `loggamma`). The real part is `log|Γ|`, which equals `std::lgamma` on the real axis. |
| `digamma(q)` | ψ = Γ′/Γ. Uses the asymptotic series after recurrence, and reflection for `Re < ½`. |
| `beta(q, b)`, `beta(b, q)` | B(q, b) = Γ(q)Γ(b)/Γ(q+b) for real `b`. Two quaternions on different axes have no well-defined beta. |
| `erf(q)`, `erfc(q)` | Maclaurin series for `|Re| < ½`. Otherwise the Laplace continued fraction for `erfc` (modified Lentz), with `erfc(−z) = 2 − erfc(z)`. |
| `riemann_zeta(q)` | ζ. Uses Borwein's η algorithm on `Re ≥ ½`, with the number of terms chosen from the precision and `|Im|`, and the functional equation on `Re < ½`. ζ(1) = ∞, ζ(0) = −½, and the trivial zeros are exact. |
| `cyl_bessel_j(n, q)` | J_n for integer `n` (an entire function). Uses the power series for `|z| ≤ 1`. Otherwise it uses Miller's backward recurrence, normalized by `1 = J₀ + 2ΣJ₂ₖ` or `cos z = J₀ + 2Σ(−1)ᵏJ₂ₖ`, whichever cancels less. |

**Branch cuts.** The complex kernel is only evaluated on the closed upper half-plane, since `|v| ≥ 0`. On a branch cut along the real axis, the result is the limit from above (imaginary part +0), which is what `std::complex` gives for `(x, +0)`. For example, `asin(2) = π/2 + 1.3169…i`, `lgamma(−2.5) = −0.0562… − 3πi`, and `log1p(−3) = log 2 + πi`. Any imaginary part of a real argument goes on the `i` axis, as with `sqrt` and `log`.

**Accuracy.** The kernels run in at least double precision, so `float` results are close to correctly rounded. Measured against mpmath at 2,800 random points (|Re| ≤ 40, |Im| ≤ 60), the worst relative error was:

| | double | long double (x87) |
|---|---|---|
| tgamma, lgamma, digamma | 2e-14 | 1.5e-16 (close to poles) |
| erf, erfc | 6e-15 | 5e-18 |
| riemann_zeta | 6e-14 | 3e-17 |
| cyl_bessel_j | 1.2e-13 (near real zeros) | 2.5e-17 |

Typical errors are a few ulp. Accuracy degrades in these cases:

* Gamma for large `|z|`: the relative error grows like `|lgamma(z)|·ε`.
* erf/erfc where `|z|²` is large: the error of `exp(−z²)` grows like `|z|²·ε`.
* ζ for `|Im s|` above a few hundred (the number of terms is capped).
* lgamma for `Re z < −4000`: it falls back to reflection, so the imaginary part may be off by a multiple of 2π.
* J_n for `|z|` beyond about 10⁵: the cost grows like `O(|z|)`.

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
* The special functions use a single complex kernel per function. The `std::` special functions that need non-integer orders or several arguments (for example `cyl_bessel_j(ν, x)` with real ν, Legendre polynomials, or elliptic integrals) are not provided.
* `std::numbers` is not used, because libc++ 19 gives only double-precision constants for `long double`. The library carries its own `long double` constants in `quat::detail`.
* The transcendental functions are not `constexpr`, because `<cmath>` is not `constexpr` until C++26. All arithmetic, `norm`, `conj`, `inverse`, `dot`, `cross` and integer `pow` are `constexpr`.
* No special IEEE handling of infinities and NaNs (C99 Annex G style) is done beyond `proj` and `abs`.
* In a `std::format` spec, the width applies to each component, not to the whole quaternion.

## Example

See [`examples/basic.cpp`](examples/basic.cpp), which is built as `example_basic` and run as a test.

## License

MIT. See [LICENSE](LICENSE).
