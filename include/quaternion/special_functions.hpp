// special_functions.hpp - inverse trigonometric/hyperbolic and special functions
// for quat::quaternion<T>.
// SPDX-License-Identifier: MIT
//
// Every function here is the quaternion extension of an analytic function f by
// axis reduction.  Write q = a + |v| u with u = v/|v| the unit vector part.
// Then
//
//     f(q) = Re f(a + i|v|) + u * Im f(a + i|v|),
//
// so f(q) lies in the same complex plane span{1, u} as q.  It commutes with q,
// agrees exactly with the complex function on the i axis, and satisfies
// f(conj q) = conj f(q) whenever f is real on the real axis.
//
// The complex kernel is only evaluated on the closed upper half plane (|v| >= 0).
// On a branch cut lying on the real axis, the result is the limit from above
// (imaginary part +0), which is what std::complex gives for (x, +0).  When v == 0
// any imaginary part is put on the i axis, as quat::sqrt and quat::log do.
//
// The kernels are evaluated in at least double precision, so float results are
// close to correctly rounded.
#ifndef QUATERNION_SPECIAL_FUNCTIONS_HPP
#define QUATERNION_SPECIAL_FUNCTIONS_HPP

#include "quaternion.hpp"

#include <algorithm>
#include <climits>
#include <cmath>
#include <complex>
#include <limits>
#include <type_traits>
#include <vector>

namespace quat {

namespace detail::sf {

template <class T>
using cplx = std::complex<T>;

// Working precision: at least double.
template <class T>
using work_t = std::conditional_t<std::is_same_v<T, float>, double, T>;

template <class T>
inline constexpr T pi = quat::detail::pi_v<T>;

template <class T>
[[nodiscard]] inline T eps() noexcept { return std::numeric_limits<T>::epsilon(); }

template <class T>
[[nodiscard]] inline T nan() noexcept { return std::numeric_limits<T>::quiet_NaN(); }

// Complex division by Smith's algorithm.  (libc++ 18/19 divides complex<long double>
// in double precision, so the kernels never use operator/ between complex values.)
template <class T>
[[nodiscard]] inline cplx<T> cdiv(const cplx<T>& a, const cplx<T>& b) noexcept {
    const T br = b.real(), bi = b.imag();
    if (br == T() && bi == T()) return {a.real() / br, a.imag() / br};
    if (std::abs(br) >= std::abs(bi)) {
        const T r = bi / br, d = br + bi * r;
        return {(a.real() + a.imag() * r) / d, (a.imag() - a.real() * r) / d};
    }
    const T r = br / bi, d = br * r + bi;
    return {(a.real() * r + a.imag()) / d, (a.imag() * r - a.real()) / d};
}
template <class T>
[[nodiscard]] inline cplx<T> cdiv(const T& a, const cplx<T>& b) noexcept { return cdiv(cplx<T>(a), b); }

// sin(pi x) and cos(pi x) for real x, with exact zeros at integers and half-integers.
template <class T>
[[nodiscard]] inline T sinpi(T x) noexcept {
    if (!std::isfinite(x)) return nan<T>();
    T r = std::fmod(x, T(2));  // (-2, 2)
    if (r <= T(-1)) r += T(2);
    else if (r > T(1)) r -= T(2);  // (-1, 1]
    if (r > T(0.5)) r = T(1) - r;
    else if (r < T(-0.5)) r = T(-1) - r;
    return std::sin(pi<T> * r);
}
template <class T>
[[nodiscard]] inline T cospi(T x) noexcept {
    if (!std::isfinite(x)) return nan<T>();
    T r = std::fmod(std::abs(x), T(2));  // [0, 2)
    if (r > T(1)) r = T(2) - r;          // [0, 1]
    return std::sin(pi<T> * (T(0.5) - r));
}
template <class T>
[[nodiscard]] inline cplx<T> sinpi(const cplx<T>& z) noexcept {
    const T py = pi<T> * z.imag();
    return {sinpi(z.real()) * std::cosh(py), cospi(z.real()) * std::sinh(py)};
}
template <class T>
[[nodiscard]] inline cplx<T> cospi(const cplx<T>& z) noexcept {
    const T py = pi<T> * z.imag();
    return {cospi(z.real()) * std::cosh(py), -sinpi(z.real()) * std::sinh(py)};
}

template <class T>
[[nodiscard]] inline bool is_nonpositive_integer(const cplx<T>& z) noexcept {
    return z.imag() == T() && z.real() <= T() && z.real() == std::floor(z.real());
}

// Even-index Bernoulli numbers B_2 .. B_40.
template <class T>
[[nodiscard]] inline T bernoulli_2k(int k) noexcept {
    static constexpr long double num[] = {
        1.0L, -1.0L, 1.0L, -1.0L, 5.0L, -691.0L, 7.0L, -3617.0L, 43867.0L, -174611.0L,
        854513.0L, -236364091.0L, 8553103.0L, -23749461029.0L, 8615841276005.0L,
        -7709321041217.0L, 2577687858367.0L, -26315271553053477373.0L, 2929993913841559.0L,
        -261082718496449122051.0L};
    static constexpr long double den[] = {6.0L, 30.0L, 42.0L, 30.0L, 66.0L, 2730.0L, 6.0L,
                                          510.0L, 798.0L, 330.0L, 138.0L, 2730.0L, 6.0L,
                                          870.0L, 14322.0L, 510.0L, 6.0L, 1919190.0L, 6.0L,
                                          13530.0L};
    return static_cast<T>(num[k - 1] / den[k - 1]);
}

// Shift point for the asymptotic (Stirling) series, and its number of terms.
// With |w| >= 20 and 16 terms the truncation error is below 1e-33.
template <class T>
inline constexpr T stirling_shift = T(20);
inline constexpr int stirling_terms = 16;
inline constexpr long max_recurrence_shift = 4000;

// ---- log-gamma ---------------------------------------------------------------
// Uses Stirling's series after shifting upwards with Gamma(z) = Gamma(z+n) / prod(z+k).
// The imaginary part is accumulated as a sum of principal logarithms, so the
// result is the analytic continuation of log Gamma from the positive real axis
// (the "loggamma" branch, continuous away from the cut along the negative real
// axis).  Its real part is always log|Gamma(z)|.
template <class T>
[[nodiscard]] inline cplx<T> lgamma(cplx<T> z) noexcept {
    if (is_nonpositive_integer(z)) return {std::numeric_limits<T>::infinity(), T()};
    if (std::isnan(z.real()) || std::isnan(z.imag())) return {nan<T>(), nan<T>()};
    cplx<T> shift_sum{};
    if (z.real() < stirling_shift<T>) {
        const T need = std::ceil(stirling_shift<T> - z.real());
        if (need > T(max_recurrence_shift)) {
            // Far left half-plane: reflection (the imaginary part may differ from
            // the continuous branch by a multiple of 2*pi).
            return std::log(pi<T>) - std::log(sinpi(z)) - lgamma(cplx<T>(T(1)) - z);
        }
        const long n = static_cast<long>(need);
        for (long k = 0; k < n; ++k) shift_sum += std::log(z + T(k));
        z += T(n);
    }
    const cplx<T> inv = cdiv(T(1), z), inv2 = inv * inv;
    cplx<T> series{}, p = inv;
    for (int k = 1; k <= stirling_terms; ++k) {
        series += bernoulli_2k<T>(k) / T(2 * k * (2 * k - 1)) * p;
        p *= inv2;
    }
    const T half_log_2pi = T(0.5) * std::log(T(2) * pi<T>);
    return (z - T(0.5)) * std::log(z) - z + half_log_2pi + series - shift_sum;
}

// ---- gamma -------------------------------------------------------------------
template <class T>
[[nodiscard]] inline cplx<T> tgamma(const cplx<T>& z) noexcept {
    if (is_nonpositive_integer(z)) {
        if (z.real() == T()) return {std::copysign(std::numeric_limits<T>::infinity(), z.real()), T()};
        return {nan<T>(), T()};
    }
    if (z.real() < T(0.5)) {
        // Reflection: Gamma(z) Gamma(1-z) = pi / sin(pi z).
        return cdiv(pi<T>, sinpi(z) * tgamma(cplx<T>(T(1)) - z));
    }
    return std::exp(lgamma(z));
}

// ---- digamma -----------------------------------------------------------------
template <class T>
[[nodiscard]] inline cplx<T> digamma(cplx<T> z) noexcept {
    if (is_nonpositive_integer(z)) return {nan<T>(), T()};
    if (z.real() < T(0.5)) {
        // psi(1-z) - psi(z) = pi cot(pi z)
        return digamma(cplx<T>(T(1)) - z) - cdiv(pi<T> * cospi(z), sinpi(z));
    }
    cplx<T> shift_sum{};
    while (z.real() < stirling_shift<T>) {
        shift_sum += cdiv(T(1), z);
        z += T(1);
    }
    const cplx<T> inv = cdiv(T(1), z), inv2 = inv * inv;
    cplx<T> series{}, p = inv2;
    for (int k = 1; k <= stirling_terms; ++k) {
        series += bernoulli_2k<T>(k) / T(2 * k) * p;
        p *= inv2;
    }
    return std::log(z) - T(0.5) * inv - series - shift_sum;
}

// ---- elementary helpers --------------------------------------------------------
template <class T>
[[nodiscard]] inline cplx<T> expm1(const cplx<T>& z) noexcept {
    const T x = z.real(), y = z.imag();
    if (y == T()) return {std::expm1(x), y};
    const T s = std::sin(y / T(2));
    return {std::expm1(x) * std::cos(y) - T(2) * s * s, std::exp(x) * std::sin(y)};
}
template <class T>
[[nodiscard]] inline cplx<T> log1p(const cplx<T>& z) noexcept {
    const T x = z.real(), y = z.imag();
    if (y == T() && x >= T(-1)) return {std::log1p(x), y};
    const T ax = std::abs(x), ay = std::abs(y);
    if (std::max(ax, ay) > T(0.5)) return std::log(T(1) + z);
    // |1+z|^2 - 1 = x(2+x) + y^2 without cancellation for small z.
    return {T(0.5) * std::log1p(x * (T(2) + x) + y * y), std::atan2(y, T(1) + x)};
}
template <class T>
[[nodiscard]] inline cplx<T> cbrt(const cplx<T>& z) noexcept {
    if (z.imag() == T() && z.real() >= T()) return {std::cbrt(z.real()), z.imag()};
    return std::polar(std::cbrt(std::abs(z)), std::arg(z) / T(3));
}

// ---- error functions -------------------------------------------------------------
// Maclaurin series of erf.  Used for |Re z| < 1/2, where the relative cancellation
// is bounded by about exp(2 Re(z)^2) < e.
template <class T>
[[nodiscard]] inline cplx<T> erf_series(const cplx<T>& z) noexcept {
    const cplx<T> z2 = z * z;
    cplx<T> term = z, sum = z;
    for (int n = 1; n < 200000; ++n) {
        term *= -z2 / T(n);
        const cplx<T> t = term / T(2 * n + 1);
        sum += t;
        if (std::abs(t) <= eps<T>() * std::abs(sum) / 4) break;
    }
    return sum * (T(2) / std::sqrt(pi<T>));
}
// Continued fraction (Laplace/Legendre) for Re z > 0:
//   erfc z = exp(-z^2)/sqrt(pi) * 1/(z + (1/2)/(z + 1/(z + (3/2)/(z + ...))))
// evaluated by the modified Lentz algorithm.
template <class T>
[[nodiscard]] inline cplx<T> erfc_cf(const cplx<T>& z) noexcept {
    const T tiny = std::numeric_limits<T>::min() * T(1e10);
    cplx<T> f = z, C = z, D{};
    if (std::abs(f) < tiny) f = C = tiny;
    for (int k = 1; k < 200000; ++k) {
        const T a = T(k) / T(2);
        D = z + a * D;
        if (std::abs(D) < tiny) D = tiny;
        D = cdiv(T(1), D);
        C = z + cdiv(a, C);
        if (std::abs(C) < tiny) C = tiny;
        const cplx<T> delta = C * D;
        f *= delta;
        if (std::abs(delta - T(1)) <= eps<T>() / 2) break;
    }
    return cdiv(std::exp(-z * z), std::sqrt(pi<T>) * f);
}
template <class T>
[[nodiscard]] inline cplx<T> erf(const cplx<T>& z) noexcept {
    if (std::isnan(z.real()) || std::isnan(z.imag())) return {nan<T>(), nan<T>()};
    if (z.imag() == T() && std::isinf(z.real())) return {std::copysign(T(1), z.real()), z.imag()};
    if (std::abs(z.real()) < T(0.5)) return erf_series(z);
    if (z.real() > T()) return T(1) - erfc_cf(z);
    return erfc_cf(-z) - T(1);
}
template <class T>
[[nodiscard]] inline cplx<T> erfc(const cplx<T>& z) noexcept {
    if (std::isnan(z.real()) || std::isnan(z.imag())) return {nan<T>(), nan<T>()};
    if (z.imag() == T() && std::isinf(z.real())) return {z.real() > 0 ? T() : T(2), z.imag()};
    if (z.real() >= T(0.5)) return erfc_cf(z);
    if (z.real() <= T(-0.5)) return T(2) - erfc_cf(-z);
    return T(1) - erf_series(z);
}

// ---- Riemann zeta ------------------------------------------------------------------
// Borwein's algorithm 2 for the Dirichlet eta function on Re s >= 1/2:
// zeta(s) = eta(s) / (1 - 2^(1-s)).  The number of terms is chosen from the target
// precision and |Im s|, so accuracy holds up to |Im s| of several hundred.
// Re s < 1/2 uses the functional equation.
template <class T>
[[nodiscard]] inline cplx<T> riemann_zeta(const cplx<T>& s) noexcept {
    if (std::isnan(s.real()) || std::isnan(s.imag())) return {nan<T>(), nan<T>()};
    if (s.imag() == T()) {
        if (s.real() == T(1)) return {std::numeric_limits<T>::infinity(), T()};
        if (s.real() == T()) return {T(-0.5), T()};
        if (s.real() < T() && std::fmod(s.real(), T(2)) == T()) return {T(), T()};  // trivial zeros
        if (s.real() == std::numeric_limits<T>::infinity()) return {T(1), T()};
    }
    if (s.real() < T(0.5)) {
        // zeta(s) = 2^s pi^(s-1) sin(pi s/2) Gamma(1-s) zeta(1-s), in log form to avoid overflow.
        const cplx<T> one_minus_s = cplx<T>(T(1)) - s;
        const cplx<T> lg = s * std::log(T(2)) + (s - T(1)) * std::log(pi<T>) + lgamma(one_minus_s);
        return std::exp(lg) * sinpi(s / T(2)) * riemann_zeta(one_minus_s);
    }
    const T t = std::abs(s.imag());
    const T target = T(std::numeric_limits<T>::digits) * std::log(T(2)) + T(3);
    const T rate = std::log(T(3) + std::sqrt(T(8)));
    const long n = std::clamp(static_cast<long>(std::ceil(
                                  (target + pi<T> * t / T(2) + std::log1p(T(2) * t)) / rate)) + 2,
                              8L, 6000L);
    // d_k = n * sum_{i<=k} (n+i-1)! 4^i / ((n-i)! (2i)!), computed in long double
    // relative to d_n through the recursion a_i / a_{i-1}.
    using L = long double;
    cplx<T> eta{};
    {
        // a_i values grow like (3+sqrt 8)^i; normalise by d_n at the end.
        std::vector<L> d(static_cast<std::size_t>(n) + 1);
        L a = L(1) / L(n), acc = a;
        d[0] = L(n) * acc;
        for (long i = 1; i <= n; ++i) {
            a *= L(4) * L(n + i - 1) * L(n - i + 1) / (L(2 * i) * L(2 * i - 1));
            acc += a;
            d[static_cast<std::size_t>(i)] = L(n) * acc;
        }
        const L dn = d[static_cast<std::size_t>(n)];
        for (long k = 0; k < n; ++k) {
            const T coeff = static_cast<T>((dn - d[static_cast<std::size_t>(k)]) / dn);
            const cplx<T> term = coeff * std::exp(-s * std::log(T(k + 1)));
            eta += (k % 2 == 0) ? term : -term;
        }
    }
    // 1 - 2^(1-s) = -expm1((1-s) ln 2), accurate near s = 1.
    const cplx<T> denom = -expm1((cplx<T>(T(1)) - s) * std::log(T(2)));
    return cdiv(eta, denom);
}

// ---- Bessel J_n, integer order ----------------------------------------------------
template <class T>
[[nodiscard]] inline cplx<T> cyl_bessel_j(int n, const cplx<T>& z) noexcept {
    if (std::isnan(z.real()) || std::isnan(z.imag())) return {nan<T>(), nan<T>()};
    if (n < 0) {
        const cplx<T> r = cyl_bessel_j(n == INT_MIN ? INT_MAX : -n, z);  // J_{-n} = (-1)^n J_n
        return (n % 2 == 0) ? r : -r;
    }
    if (z == cplx<T>()) return {n == 0 ? T(1) : T(), T()};
    const T az = std::abs(z);
    if (az <= T(1)) {
        // Power series: no significant cancellation for |z| <= 1.
        const cplx<T> h = z / T(2), mh2 = -h * h;
        cplx<T> lead(T(1));
        for (int j = 1; j <= n; ++j) lead *= h / T(j);
        cplx<T> term(T(1)), sum(T(1));
        for (int k = 1; k < 1000; ++k) {
            term *= mh2 / (T(k) * T(n + k));
            sum += term;
            if (std::abs(term) <= eps<T>() * std::abs(sum) / 4) break;
        }
        return lead * sum;
    }
    // Miller's backward recurrence f_{k-1} = (2k/z) f_k - f_{k+1}, normalised by
    // either 1 = J0 + 2 sum J_{2k} or cos z = J0 + 2 sum (-1)^k J_{2k}, whichever
    // suffers less cancellation.
    const T digits = T(std::numeric_limits<T>::digits10);
    const long M0 = static_cast<long>(std::max(T(n), az) + T(10) + T(3) * digits +
                                      T(6) * std::cbrt(az));
    const long M = M0 + (M0 % 2);  // even start
    const cplx<T> two_over_z = cdiv(T(2), z);
    cplx<T> fk1{}, fk(std::numeric_limits<T>::min() * T(1e20));  // f_{M+1}, f_M
    cplx<T> s1{}, s2{}, fn{};
    T abssum{};
    const T big = T(1e100), small = T(1e-100);
    for (long k = M; k >= 1; --k) {
        const cplx<T> fkm1 = T(k) * two_over_z * fk - fk1;  // f_{k-1}
        fk1 = fk;
        fk = fkm1;
        const long idx = k - 1;
        if (idx == n) fn = fk;
        if (idx > 0 && idx % 2 == 0) {
            s1 += T(2) * fk;
            s2 += ((idx / 2) % 2 == 0 ? T(2) : T(-2)) * fk;
            abssum += T(2) * std::abs(fk);
        }
        if (std::abs(fk) > big) {
            fk *= small; fk1 *= small; s1 *= small; s2 *= small; fn *= small; abssum *= small;
        }
    }
    s1 += fk;  // f_0
    s2 += fk;
    abssum += std::abs(fk);
    if (n == 0) fn = fk;
    const T c1 = abssum / std::abs(s1), c2 = abssum / std::abs(s2);
    if (!(c2 < c1)) return cdiv(fn, s1);
    return cdiv(fn * std::cos(z), s2);
}

// std::acos with an exact pi/2 at zero (libc++ 19 returns a double-precision pi/2
// for complex<long double> there).
template <class T>
[[nodiscard]] inline cplx<T> acos(const cplx<T>& z) noexcept {
    if (z.real() == T() && z.imag() == T()) return {pi<T> / T(2), -z.imag()};
    return std::acos(z);
}

// atanh and atan with the C99/std::complex branch cuts, written in terms of log1p
// and cdiv (libc++'s complex<long double> versions lose precision through division).
//   atanh z = (log1p z - log1p(-z)) / 2                       for |z| <= 1
//   atanh z = atanh(1/z) + i copysign(pi/2, Im z)              for |z| >  1
//   atan  z = -i atanh(i z)
template <class T>
[[nodiscard]] inline cplx<T> atanh(const cplx<T>& z) noexcept {
    if (std::isnan(z.real()) || std::isnan(z.imag())) return std::atanh(z);
    if (std::abs(z) <= T(1)) return (log1p(z) - log1p(-z)) / T(2);
    const cplx<T> w = sf::atanh(cdiv(T(1), z));
    return {w.real(), w.imag() + std::copysign(pi<T> / T(2), z.imag())};
}
template <class T>
[[nodiscard]] inline cplx<T> atan(const cplx<T>& z) noexcept {
    const cplx<T> w = sf::atanh(cplx<T>(-z.imag(), z.real()));
    return {w.imag(), -w.real()};
}

// ---- the axis-reduction lift ----------------------------------------------------------
// F maps std::complex<W> -> std::complex<W>.  If real_on_real is true, the function is
// known to be real on the real axis, and a purely real argument yields a real result.
template <class T, class F>
[[nodiscard]] inline quaternion<T> lift(const quaternion<T>& q, F&& f, bool real_on_real) {
    using W = work_t<T>;
    const W v = static_cast<W>(std::hypot(q.x(), q.y(), q.z()));
    const cplx<W> r = f(cplx<W>(static_cast<W>(q.w()), v));
    if (v == W()) {
        if (real_on_real) return {static_cast<T>(r.real()), q.x(), q.y(), q.z()};
        return {static_cast<T>(r.real()), static_cast<T>(r.imag()), T(), T()};
    }
    const W s = r.imag() / v;
    return {static_cast<T>(r.real()), static_cast<T>(static_cast<W>(q.x()) * s),
            static_cast<T>(static_cast<W>(q.y()) * s), static_cast<T>(static_cast<W>(q.z()) * s)};
}

}  // namespace detail::sf

#define QUAT_SF_LIFT(NAME, EXPR, REAL)                                                   \
    template <class T>                                                                   \
    [[nodiscard]] inline quaternion<T> NAME(const quaternion<T>& q) noexcept {           \
        return detail::sf::lift(                                                         \
            q, [](const auto& z) { return EXPR; }, REAL);                                \
    }

// Inverse trigonometric and hyperbolic functions (std::complex branch cuts).
QUAT_SF_LIFT(asin, std::asin(z), false)
QUAT_SF_LIFT(acos, detail::sf::acos(z), false)
QUAT_SF_LIFT(atan, detail::sf::atan(z), false)
QUAT_SF_LIFT(asinh, std::asinh(z), false)
QUAT_SF_LIFT(acosh, std::acosh(z), false)
QUAT_SF_LIFT(atanh, detail::sf::atanh(z), false)

// Elementary functions with good accuracy near zero.
QUAT_SF_LIFT(expm1, detail::sf::expm1(z), true)
QUAT_SF_LIFT(log1p, detail::sf::log1p(z), false)
QUAT_SF_LIFT(cbrt, detail::sf::cbrt(z), false)

// Gamma family.
QUAT_SF_LIFT(tgamma, detail::sf::tgamma(z), true)
QUAT_SF_LIFT(lgamma, detail::sf::lgamma(z), false)
QUAT_SF_LIFT(digamma, detail::sf::digamma(z), true)

// Error functions.
QUAT_SF_LIFT(erf, detail::sf::erf(z), true)
QUAT_SF_LIFT(erfc, detail::sf::erfc(z), true)

// Riemann zeta function.
QUAT_SF_LIFT(riemann_zeta, detail::sf::riemann_zeta(z), true)
#undef QUAT_SF_LIFT

/// Cylinder Bessel function of the first kind J_n(q) for integer order n (entire in q).
template <class T>
[[nodiscard]] inline quaternion<T> cyl_bessel_j(int n, const quaternion<T>& q) noexcept {
    return detail::sf::lift(q, [n](const auto& z) { return detail::sf::cyl_bessel_j(n, z); }, true);
}

/// Euler beta function B(q, b) = Gamma(q) Gamma(b) / Gamma(q + b) for real b.
/// (For two quaternions on different axes the product is not well defined.)
template <class T>
[[nodiscard]] inline quaternion<T> beta(const quaternion<T>& q, const T& b) noexcept {
    using W = detail::sf::work_t<T>;
    const W wb = static_cast<W>(b);
    return detail::sf::lift(
        q,
        [wb](const std::complex<W>& z) {
            using namespace detail::sf;
            const cplx<W> cb(wb);
            if (is_nonpositive_integer(z) || is_nonpositive_integer(cb))
                return cplx<W>(nan<W>(), W());
            return std::exp(lgamma(z) + lgamma(cb) - lgamma(z + cb));
        },
        true);
}
template <class T>
[[nodiscard]] inline quaternion<T> beta(const T& a, const quaternion<T>& q) noexcept {
    return beta(q, a);
}

}  // namespace quat

#endif  // QUATERNION_SPECIAL_FUNCTIONS_HPP
