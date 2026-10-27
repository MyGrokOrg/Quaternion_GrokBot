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
#include <concepts>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <vector>

#if defined(QUATERNION_MSVC_CONTRACT_OFF)
#pragma float_control(push)  // see quaternion.hpp
#pragma fp_contract(off)
#endif

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
    return detail::lm::sin(pi<T> * r);
}
template <class T>
[[nodiscard]] inline T cospi(T x) noexcept {
    if (!std::isfinite(x)) return nan<T>();
    T r = std::fmod(std::abs(x), T(2));  // [0, 2)
    if (r > T(1)) r = T(2) - r;          // [0, 1]
    return detail::lm::sin(pi<T> * (T(0.5) - r));
}
template <class T>
[[nodiscard]] inline cplx<T> sinpi(const cplx<T>& z) noexcept {
    const T py = pi<T> * z.imag();
    return {sinpi(z.real()) * detail::lm::cosh(py), cospi(z.real()) * detail::lm::sinh(py)};
}
template <class T>
[[nodiscard]] inline cplx<T> cospi(const cplx<T>& z) noexcept {
    const T py = pi<T> * z.imag();
    return {cospi(z.real()) * detail::lm::cosh(py), -sinpi(z.real()) * detail::lm::sinh(py)};
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

// ---- double-T ("double-double") arithmetic ---------------------------------------------
// A value hi + lo with |lo| <= ulp(hi)/2, giving about twice the precision of T.  Used
// where a result is an exponential of a large argument (Gamma, erfc, zeta), so that the
// argument is carried with an absolute error far below ulp(1).  Products use std::fma.
template <class T>
struct dd {
    T hi{}, lo{};
};
template <class T>
[[nodiscard]] inline dd<T> two_sum(T a, T b) noexcept {
    const T s = a + b, bb = s - a;
    return {s, (a - (s - bb)) + (b - bb)};
}
template <class T>
[[nodiscard]] inline dd<T> fast_two_sum(T a, T b) noexcept {
    const T s = a + b;
    return {s, b - (s - a)};
}
template <class T>
[[nodiscard]] inline dd<T> two_prod(T a, T b) noexcept {
    const T p = a * b;
    if constexpr (std::numeric_limits<T>::digits == 64) {
        // x87 extended: std::fma is a slow software routine, so use Dekker's exact product
        // (Veltkamp split at 2^32 + 1); exact unless a or b is within 2^32 of overflow.
        // (2^16000 through pow2: MSVC rejects the literal 0x1p+16000L even in a discarded
        // branch, since its long double is double.)
        constexpr T big = detail::cx::pow2<T>(16000);
        if (std::abs(a) < big && std::abs(b) < big) {
            constexpr T split = T(4294967297.0L);
            const T ta = split * a, ah = ta - (ta - a), al = a - ah;
            const T tb = split * b, bh = tb - (tb - b), bl = b - bh;
            return {p, ((ah * bh - p) + ah * bl + al * bh) + al * bl};
        }
    }
    return {p, std::fma(a, b, -p)};
}
template <class T>
[[nodiscard]] inline dd<T> operator-(const dd<T>& a) noexcept { return {-a.hi, -a.lo}; }
template <class T>
[[nodiscard]] inline dd<T> operator+(const dd<T>& a, const dd<T>& b) noexcept {
    dd<T> s = two_sum(a.hi, b.hi);
    s.lo += a.lo + b.lo;
    return fast_two_sum(s.hi, s.lo);
}
template <class T>
[[nodiscard]] inline dd<T> operator+(const dd<T>& a, T b) noexcept {
    dd<T> s = two_sum(a.hi, b);
    s.lo += a.lo;
    return fast_two_sum(s.hi, s.lo);
}
template <class T>
[[nodiscard]] inline dd<T> operator-(const dd<T>& a, const dd<T>& b) noexcept { return a + (-b); }
template <class T>
[[nodiscard]] inline dd<T> operator-(const dd<T>& a, T b) noexcept { return a + (-b); }
template <class T>
[[nodiscard]] inline dd<T> operator*(const dd<T>& a, const dd<T>& b) noexcept {
    dd<T> p = two_prod(a.hi, b.hi);
    p.lo += a.hi * b.lo + a.lo * b.hi;
    return fast_two_sum(p.hi, p.lo);
}
template <class T>
[[nodiscard]] inline dd<T> operator*(const dd<T>& a, T b) noexcept {
    dd<T> p = two_prod(a.hi, b);
    p.lo += a.lo * b;
    return fast_two_sum(p.hi, p.lo);
}
template <class T>
[[nodiscard]] inline dd<T> operator/(const dd<T>& a, const dd<T>& b) noexcept {
    const T q1 = a.hi / b.hi;
    dd<T> r = a - b * q1;
    const T q2 = r.hi / b.hi;
    r = r - b * q2;
    const T q3 = r.hi / b.hi;
    return fast_two_sum(q1, q2) + q3;
}
template <class T>
[[nodiscard]] inline dd<T> ldexp(const dd<T>& a, int e) noexcept {
    return {c_scalbn(a.hi, e), c_scalbn(a.lo, e)};  // = std::ldexp (exact helper on Apple x86_64)
}

enum class ddc { ln2, pi, half_log_2pi, log_pi, two_over_sqrt_pi };
// Constants split as hi + lo for 64-bit (x87) and 53-bit significands.
template <class T>
[[nodiscard]] inline dd<T> dd_const(ddc c) noexcept {
    if constexpr (std::numeric_limits<T>::digits == 64) {
        switch (c) {
            case ddc::ln2: return {T(0x2c5c85fdf473de6bp-62L), T(-0xd871319ff0342543p-130L)};
            case ddc::pi: return {T(0xc90fdaa22168c235p-62L), T(-0xece675d1fc8f8cbbp-128L)};
            case ddc::half_log_2pi: return {T(0xeb3f8e4325f5a535p-64L), T(-0x1ad0dbffaef9b7f7p-126L)};
            case ddc::log_pi: return {T(0x928682473d0de85fp-63L), T(-0x140d5272af7816cdp-125L)};
            case ddc::two_over_sqrt_pi: return {T(0x906eba8214db688dp-63L), T(0x71d48a7f6bfec344p-127L)};
        }
    } else if constexpr (std::numeric_limits<T>::digits == 53) {
        switch (c) {
            case ddc::ln2: return {T(0x1.62e42fefa39efp-1), T(0x1.abc9e3b39803fp-56)};
            case ddc::pi: return {T(0x1.921fb54442d18p+1), T(0x1.1a62633145c07p-53)};
            case ddc::half_log_2pi: return {T(0x1.d67f1c864beb5p-1), T(-0x1.65b5a1b7ff5dfp-55)};
            case ddc::log_pi: return {T(0x1.250d048e7a1bdp+0), T(0x1.7abf2ad8d5088p-57)};
            case ddc::two_over_sqrt_pi: return {T(0x1.20dd750429b6dp+0), T(0x1.1ae3a914fed80p-56)};
        }
    }
    // Other formats: no extra precision.
    switch (c) {
        case ddc::ln2: return {std::log(T(2)), T()};
        case ddc::pi: return {pi<T>, T()};
        case ddc::half_log_2pi: return {T(0.5) * std::log(T(2) * pi<T>), T()};
        case ddc::two_over_sqrt_pi: return {T(2) / detail::lm::sqrt(pi<T>), T()};
        default: return {std::log(pi<T>), T()};
    }
}

// exp of a double-T value: 2^k reduction, a further 2^-8 scaling, Taylor series for
// expm1, then eight squarings of 1 + expm1.
template <class T>
[[nodiscard]] inline dd<T> exp_dd(const dd<T>& x) noexcept {
    if (std::isnan(x.hi)) return {x.hi, T()};
    if (x.hi > std::log(std::numeric_limits<T>::max())) return {std::numeric_limits<T>::infinity(), T()};
    if (x.hi < std::log(std::numeric_limits<T>::denorm_min()) - T(1)) return {};
    const dd<T> ln2 = dd_const<T>(ddc::ln2);
    const T k = std::nearbyint(x.hi / ln2.hi);
    const dd<T> r = ldexp(x - ln2 * k, -8);
    dd<T> s = r, t = r;
    for (int i = 2; i <= 14; ++i) {
        t = t * r / dd<T>{T(i), T()};
        s = s + t;
    }
    for (int i = 0; i < 8; ++i) s = s * (s + T(2));
    return ldexp(s + T(1), static_cast<int>(k));
}
// log of a positive double-T value by Newton's method on exp_dd.
template <class T>
[[nodiscard]] inline dd<T> log_dd(const dd<T>& x) noexcept {
    dd<T> y{std::log(x.hi), T()};
    for (int i = 0; i < 2; ++i) y = y + x * exp_dd(-y) - T(1);
    return y;
}
// sin and cos of a real T value with |a| <= 4, to double-T precision.
template <class T>
inline void sincos_dd(T a, dd<T>& s, dd<T>& c) noexcept {
    const dd<T> pio2 = ldexp(dd_const<T>(ddc::pi), -1);
    const T j = std::nearbyint(a / pio2.hi);
    const dd<T> r = dd<T>{a, T()} - pio2 * j;
    const dd<T> r2 = r * r;
    dd<T> ts = r, tc{T(1), T()};
    s = ts;
    c = tc;
    for (int i = 1; i <= 16; ++i) {
        ts = -(ts * r2) / dd<T>{T((2 * i) * (2 * i + 1)), T()};
        tc = -(tc * r2) / dd<T>{T((2 * i - 1) * (2 * i)), T()};
        s = s + ts;
        c = c + tc;
    }
    const int q = static_cast<int>(j) & 3;
    const dd<T> s0 = s, c0 = c;
    if (q == 1) s = c0, c = -s0;
    else if (q == 2) s = -s0, c = -c0;
    else if (q == 3) s = -c0, c = s0;
}

template <class T>
struct cdd {
    dd<T> re, im;
};
template <class T>
[[nodiscard]] inline cdd<T> operator+(const cdd<T>& a, const cdd<T>& b) noexcept {
    return {a.re + b.re, a.im + b.im};
}
template <class T>
[[nodiscard]] inline cdd<T> operator-(const cdd<T>& a, const cdd<T>& b) noexcept {
    return {a.re - b.re, a.im - b.im};
}
template <class T>
[[nodiscard]] inline cdd<T> operator*(const cdd<T>& a, const cdd<T>& b) noexcept {
    return {a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re};
}
template <class T>
[[nodiscard]] inline cdd<T> operator*(const cdd<T>& a, const dd<T>& b) noexcept {
    return {a.re * b, a.im * b};
}
template <class T>
[[nodiscard]] inline cdd<T> to_cdd(const cplx<T>& z) noexcept {
    return {{z.real(), T()}, {z.imag(), T()}};
}
template <class T>
[[nodiscard]] inline cplx<T> to_cplx(const cdd<T>& z) noexcept {
    return {z.re.hi, z.im.hi};
}
// exp(hi + lo) = exp(hi) (1 + lo) for a double-T exponent.
template <class T>
[[nodiscard]] inline cplx<T> exp_split(const cdd<T>& e) noexcept {
    const cplx<T> h = detail::lm::exp(cplx<T>(e.re.hi, e.im.hi));
    const T lr = e.re.lo, li = e.im.lo;
    return {h.real() * (T(1) + lr) - h.imag() * li, h.imag() * (T(1) + lr) + h.real() * li};
}
// Principal log of x + i y (x, y double-T) to double-T precision.
template <class T>
[[nodiscard]] inline cdd<T> clog_dd(const dd<T>& x, const dd<T>& y) noexcept {
    const T m = std::max(std::abs(x.hi), std::abs(y.hi));
    const int e = std::ilogb(m);
    const dd<T> xs = ldexp(x, -e), ys = ldexp(y, -e);
    const dd<T> r2 = xs * xs + ys * ys;
    const dd<T> re = ldexp(log_dd(r2), -1) + dd_const<T>(ddc::ln2) * T(e);
    const T t0 = std::atan2(y.hi, x.hi);
    dd<T> s, c;
    sincos_dd(t0, s, c);
    // tan(theta - t0) = (y cos t0 - x sin t0) / (x cos t0 + y sin t0)
    const dd<T> d = (ys * c - xs * s) / (xs * c + ys * s);
    return {re, dd<T>{t0, T()} + d};
}

// Plain T evaluation, used only for |z| > 1e30 where double-T products could overflow.
template <class T>
[[nodiscard]] inline cplx<T> lgamma_plain(cplx<T> z) noexcept {
    if (is_nonpositive_integer(z)) return {std::numeric_limits<T>::infinity(), T()};
    if (std::isnan(z.real()) || std::isnan(z.imag())) return {nan<T>(), nan<T>()};
    cplx<T> shift_sum{};
    if (z.real() < stirling_shift<T>) {
        const T need = std::ceil(stirling_shift<T> - z.real());
        if (need > T(max_recurrence_shift)) {
            // Far left half-plane: reflection (the imaginary part may differ from
            // the continuous branch by a multiple of 2*pi).
            return std::log(pi<T>) - std::log(sinpi(z)) - lgamma_plain(cplx<T>(T(1)) - z);
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

// ---- log-gamma ---------------------------------------------------------------
// Uses Stirling's series after shifting upwards with Gamma(z) = Gamma(z+n) / prod(z+k).
// The imaginary part is accumulated as a sum of principal logarithms, so the
// result is the analytic continuation of log Gamma from the positive real axis
// (the "loggamma" branch, continuous away from the cut along the negative real
// axis).  Its real part is always log|Gamma(z)|.
//
// lgamma_dd evaluates this in double-T arithmetic (the shifted argument, the logs and
// (z - 1/2) log z - z are carried as hi + lo; only the small Stirling correction is in
// T), so exp(lgamma) loses no accuracy to the size of log Gamma.  The real part of the
// argument is double-T as well, so that 1 - z is exact in the reflection formula.
template <class T>
[[nodiscard]] inline cdd<T> lgamma_dd(dd<T> x, T y) noexcept {
    const dd<T> yd{y, T()};
    cdd<T> shift_sum{};
    if (x.hi < stirling_shift<T>) {
        const long n = static_cast<long>(std::ceil(stirling_shift<T> - x.hi));
        for (long k = 0; k < n; ++k) shift_sum = shift_sum + clog_dd(x + T(k), yd);
        x = x + T(n);
    }
    const cplx<T> z(x.hi, y);
    const cplx<T> inv = cdiv(T(1), z), inv2 = inv * inv;
    cplx<T> series{}, p = inv;
    for (int k = 1; k <= stirling_terms; ++k) {
        series += bernoulli_2k<T>(k) / T(2 * k * (2 * k - 1)) * p;
        p *= inv2;
    }
    const cdd<T> lz = clog_dd(x, yd);
    const cdd<T> zmh{x - T(0.5), yd};
    cdd<T> r = zmh * lz - cdd<T>{x, yd};
    r.re = r.re + dd_const<T>(ddc::half_log_2pi);
    r = r + to_cdd(series) - shift_sum;
    return r;
}

// sin(pi z) = exp(E) * P with E real double-T, so that sin(pi z) can be used in log
// form for large |Im z| without overflow:
//   sin(pi z) = e^{pi |y|} e^{-+ i pi x} (w - 1) (-+ i/2),  w = e^{+-2 pi i x} e^{-2 pi |y|}.
template <class T>
struct sinpi_exp {
    dd<T> E;
    cplx<T> P;
};
template <class T>
[[nodiscard]] inline sinpi_exp<T> sinpi_split(const cplx<T>& z) noexcept {
    const T x = z.real(), y = z.imag();
    if (std::abs(y) <= T(1)) return {{}, sinpi(z)};
    const T ay = std::abs(y);
    const T ew = detail::lm::exp(T(-2) * pi<T> * ay);
    const cplx<T> w(ew * cospi(T(2) * x), ew * sinpi(T(2) * x));
    const cplx<T> ph(cospi(x), -sinpi(x));  // e^{-i pi x}
    cplx<T> P = ph * (w - T(1)) * cplx<T>(T(), T(-0.5));
    if (y < T()) P = std::conj(P);
    return {dd_const<T>(ddc::pi) * ay, P};
}

// Reflection for Re z far below zero, where the upward shift would be too long:
//   log Gamma(z) = log pi - log sin(pi z) - log Gamma(1 - z) + 2 pi i m.
// The integer m restoring the continuous branch is found from Stirling's formula
// S(z) = (z - 1/2) log z - z + log(2 pi)/2, whose imaginary part differs from the true
// one by the imaginary part of -log(1 - e^{2 pi i z}), which is below pi/2 in size.
template <class T>
[[nodiscard]] inline cplx<T> lgamma_reflect(cplx<T> z) noexcept {
    const bool lower = z.imag() < T() || (z.imag() == T() && std::signbit(z.imag()));
    if (lower) z = {z.real(), -z.imag()};
    const sinpi_exp<T> sp = sinpi_split(z);
    const cdd<T> lg1 = lgamma_dd(two_sum(T(1), -z.real()), -z.imag());
    const cplx<T> ls = sp.E.hi + std::log(sp.P);
    const cplx<T> R = cplx<T>(std::log(pi<T>) - lg1.re.hi, -lg1.im.hi) - ls;
    const cplx<T> S = (z - T(0.5)) * std::log(z) - z;
    const T m = std::nearbyint((S.imag() - R.imag()) / (T(2) * pi<T>));
    const dd<T> im = dd<T>{-lg1.im.hi, -lg1.im.lo} - dd<T>{ls.imag(), T()} +
                     dd_const<T>(ddc::pi) * (T(2) * m);
    const cplx<T> res(R.real(), im.hi);
    return lower ? std::conj(res) : res;
}

template <class T>
[[nodiscard]] inline cplx<T> lgamma(const cplx<T>& z) noexcept {
    if (is_nonpositive_integer(z)) return {std::numeric_limits<T>::infinity(), T()};
    if (std::isnan(z.real()) || std::isnan(z.imag())) return {nan<T>(), nan<T>()};
    if (std::abs(z) > T(1e30)) return lgamma_plain(z);
    if (z.real() < stirling_shift<T> &&
        std::ceil(stirling_shift<T> - z.real()) > T(max_recurrence_shift))
        return lgamma_reflect(z);
    return to_cplx(lgamma_dd(dd<T>{z.real(), T()}, z.imag()));
}

// ---- gamma -------------------------------------------------------------------
// Gamma(z) = exp(log Gamma(z)) with log Gamma in double-T, so the only error growth with
// |z| is the inherent conditioning |z psi(z)| eps of the argument itself.  Re z < 1/2
// uses Gamma(z) = pi / (sin(pi z) Gamma(1 - z)) with sin(pi z) in split exponential form.
template <class T>
[[nodiscard]] inline cplx<T> tgamma(const cplx<T>& z) noexcept {
    if (is_nonpositive_integer(z)) {
        if (z.real() == T()) return {std::copysign(std::numeric_limits<T>::infinity(), z.real()), T()};
        return {nan<T>(), T()};
    }
    if (std::isnan(z.real()) || std::isnan(z.imag())) return {nan<T>(), nan<T>()};
    if (std::abs(z) > T(1e30)) return detail::lm::exp(lgamma(z));
    if (z.real() >= T(0.5)) return exp_split(lgamma_dd(dd<T>{z.real(), T()}, z.imag()));
    const sinpi_exp<T> sp = sinpi_split(z);
    const cdd<T> lg1 = lgamma_dd(two_sum(T(1), -z.real()), -z.imag());
    cdd<T> e{dd_const<T>(ddc::log_pi) - lg1.re - sp.E, -lg1.im};
    return cdiv(exp_split(e), sp.P);
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
    if (y == T()) return {detail::lm::expm1(x), y};
    const T s = detail::lm::sin(y / T(2));
    return {detail::lm::expm1(x) * detail::lm::cos(y) - T(2) * s * s, detail::lm::exp(x) * detail::lm::sin(y)};
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
    return detail::lm::polar(std::cbrt(std::abs(z)), std::arg(z) / T(3));
}

// ---- error functions -------------------------------------------------------------
// Maclaurin series of erf.  Used for |Re z| < 1/2, where the relative cancellation
// is bounded by about exp(2 Re(z)^2) < e.  For large |Im z| the series has ~e|z|^2
// terms, so z^2, the terms and the sum are carried in double-T; otherwise the
// rounding of z^2 and of each term would accumulate to ~|z|^2 eps.
template <class T>
[[nodiscard]] inline cdd<T> erf_series_dd(const cplx<T>& z) noexcept {
    const T x = z.real(), y = z.imag();
    const cdd<T> mw{-(two_prod(x, x) - two_prod(y, y)), -(two_prod(x, y) * T(2))};
    cdd<T> term = to_cdd(z), sum = term;
    for (int n = 1; n < 200000; ++n) {
        term = term * mw;
        const dd<T> inv = dd<T>{T(1), T()} / dd<T>{T(n), T()};
        term = term * inv;
        const dd<T> inv2 = dd<T>{T(1), T()} / dd<T>{T(2 * n + 1), T()};
        const cdd<T> t = term * inv2;
        sum = sum + t;
        if (std::hypot(t.re.hi, t.im.hi) <= eps<T>() * eps<T>() * std::hypot(sum.re.hi, sum.im.hi)) break;
    }
    return sum * dd_const<T>(ddc::two_over_sqrt_pi);
}
template <class T>
[[nodiscard]] inline cplx<T> erf_series(const cplx<T>& z) noexcept {
    return to_cplx(erf_series_dd(z));
}
// erfc = 1 - erf with erf from the double-T series, rounded once.
template <class T>
[[nodiscard]] inline cplx<T> erfc_series(const cplx<T>& z) noexcept {
    const cdd<T> e = erf_series_dd(z);
    return {(-e.re + T(1)).hi, -e.im.hi};
}
// The Laplace continued fraction needs ~700 Lentz steps near |z| = 1/2 (rounding errors
// grow with the step count, and a contracted (FMA) build roughly doubles them: up to 47.9 eps
// for erfc at 0.4 <= |Re z| <= 1.3 with gcc -march=haswell, 19.4 eps without FMA, before the
// FMA audit; 12.3 / 8.46 eps after, 300 points, seed 1, fuzz/fma_audit.py --sites).  Inside
// |Re z| < 5/4, |z| < 4 the double-T series is accurate (cancellation <= e^{2 Re(z)^2}
// eps^2) and faster.
template <class T>
[[nodiscard]] inline bool erf_use_series(const cplx<T>& z) noexcept {
    const T x = std::abs(z.real());
    return x < T(0.5) || (x < T(1.25) && x * x + z.imag() * z.imag() < T(16));
}
// exp(-z^2) with z^2 = (x^2 - y^2) + 2ixy formed exactly in double-T (fma products), so
// that a large |z|^2 does not amplify the rounding of z^2 itself.
template <class T>
[[nodiscard]] inline cplx<T> exp_neg_sq(const cplx<T>& z) noexcept {
    const T x = z.real(), y = z.imag();
    const dd<T> re = two_prod(x, x) - two_prod(y, y);
    const dd<T> im = two_prod(x, y) * T(2);
    if (!std::isfinite(re.hi) || !std::isfinite(im.hi)) return detail::lm::exp(-z * z);
    return exp_split(cdd<T>{-re, -im});
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
    return cdiv(exp_neg_sq(z), detail::lm::sqrt(pi<T>) * f);
}
template <class T>
[[nodiscard]] inline cplx<T> erf(const cplx<T>& z) noexcept {
    if (std::isnan(z.real()) || std::isnan(z.imag())) return {nan<T>(), nan<T>()};
    if (z.imag() == T() && std::isinf(z.real())) return {std::copysign(T(1), z.real()), z.imag()};
    if (erf_use_series(z)) return erf_series(z);
    if (z.real() > T()) return T(1) - erfc_cf(z);
    return erfc_cf(-z) - T(1);
}
template <class T>
[[nodiscard]] inline cplx<T> erfc(const cplx<T>& z) noexcept {
    if (std::isnan(z.real()) || std::isnan(z.imag())) return {nan<T>(), nan<T>()};
    if (z.imag() == T() && std::isinf(z.real())) return {z.real() > 0 ? T() : T(2), z.imag()};
    if (erf_use_series(z)) return erfc_series(z);
    if (z.real() > T()) return erfc_cf(z);
    return T(2) - erfc_cf(-z);
}

// Euler-Maclaurin summation for zeta(s), used for large |Im s| (no cap on |Im s|):
//   zeta(s) = sum_{k<N} k^{-s} + N^{1-s}/(s-1) + N^{-s}/2
//             + sum_j c_j s(s+1)...(s+2j-2) N^{-s-2j+1},   c_j = B_2j/(2j)!.
// For j > 20, c_j = 2 (-1)^{j+1} zeta(2j) / (2 pi)^{2j} is formed directly, so up to 80 tail terms
// are available; with N = 0.3 |s| + 30 they shrink by about (|s|/(2 pi N))^2 < 0.29 each.
// Each k^{-s} = exp(-s log k) uses log k in double-T (the phase Im(s) log k is large); the
// logarithms come from a smallest-prime-factor sieve, log k = log p + log(k/p) (exact in
// double-T), so only primes need a log evaluation.  The sum is accumulated in double-T.
template <class T>
[[nodiscard]] inline cplx<T> zeta_em(const cplx<T>& s) noexcept {
    // N below must fit a long (an infinite s once built a vector of 9e18 entries)
    if (!(std::abs(s) < T(1e15))) return {nan<T>(), nan<T>()};
    const long N = static_cast<long>(std::ceil(T(0.3) * std::abs(s))) + 30;
    const long cap = std::min<long>(N, 1L << 20);  // sieve/table size (memory bound)
    std::vector<std::int32_t> spf(static_cast<std::size_t>(cap) + 1, 0);
    std::vector<dd<T>> lg(static_cast<std::size_t>(cap) + 1);
    for (long i = 2; i <= cap; ++i) {
        if (spf[static_cast<std::size_t>(i)] == 0) {
            for (long j = i; j <= cap; j += i)
                if (spf[static_cast<std::size_t>(j)] == 0) spf[static_cast<std::size_t>(j)] = static_cast<std::int32_t>(i);
        }
    }
    auto logk = [&](long k) -> dd<T> {
        if (k > cap) return log_dd(dd<T>{T(k), T()});
        return lg[static_cast<std::size_t>(k)];
    };
    auto kpow = [&](const dd<T>& lk) { return exp_split(cdd<T>{-(lk * s.real()), -(lk * s.imag())}); };
    cdd<T> sum{{T(1), T()}, {T(), T()}};
    for (long k = 2; k <= N; ++k) {
        if (k <= cap) {
            const long p = spf[static_cast<std::size_t>(k)];
            lg[static_cast<std::size_t>(k)] =
                p == k ? log_dd(dd<T>{T(k), T()}) : lg[static_cast<std::size_t>(p)] + lg[static_cast<std::size_t>(k / p)];
        }
        if (k < N) sum = sum + to_cdd(kpow(logk(k)));
    }
    const cplx<T> nps = kpow(logk(N));  // N^{-s}
    cplx<T> tail = cdiv(T(N) * nps, s - T(1)) + nps / T(2);
    cplx<T> pp = s * (nps / T(N));  // s(s+1)...(s+2j-2) N^{-s-2j+1}, kept as one product
    const T twopi2 = T(4) * pi<T> * pi<T>;
    T pj = T(1);                        // (2 pi)^{-2j}
    T fact = T(1);                      // (2j)!
    T prev = std::numeric_limits<T>::infinity();
    for (int j = 1; j <= 80; ++j) {
        pj /= twopi2;
        fact *= T(2 * j - 1) * T(2 * j);
        T cj;
        if (j <= 20) {
            cj = bernoulli_2k<T>(j) / fact;
        } else {
            T z2j = T(1);  // zeta(2j), 2j > 40: a few terms suffice
            for (int m = 2; m < 16; ++m) {
                const T t = detail::lm::pow(T(m), T(-2 * j));
                z2j += t;
                if (t < eps<T>() * eps<T>()) break;
            }
            cj = (j % 2 ? T(2) : T(-2)) * z2j * pj;
        }
        const cplx<T> term = cj * pp;
        const T at = std::abs(term);
        tail += term;
        if (at <= eps<T>() * eps<T>() * std::abs(tail) || at > prev) break;
        prev = at;
        pp *= ((s + T(2 * j - 1)) / T(N)) * ((s + T(2 * j)) / T(N));
    }
    return to_cplx(sum + to_cdd(tail));
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
    if (!std::isfinite(s.real()) || !std::isfinite(s.imag())) return {nan<T>(), nan<T>()};  // the public function has the limits
    if (std::abs(s) < T(0.25)) {
        // Taylor series about 0.  The functional equation would need zeta(1-s) with 1-s
        // rounded, an error of eps/|s| next to the pole (2.6e6 ulp at s = -1e-7 in double).
        static constexpr long double c[] = {
        -0.5L,
        -0.91893853320467274178032973640561764L,
        -1.00317822795429242560505001336498022L,
        -1.00078519447704240796017680222772921L,
        -0.999879299500571164957800813655875236L,
        -1.00000194089632045603779988198163183L,
        -1.00000130114601395962431150487297972L,
        -0.999999831384173610779930217058015407L,
        -1.00000000576467597994939441606374166L,
        -1.00000000091101648923141657092186742L,
        -0.999999999850299240580988626479279943L,
        -1.00000000000940689566566617690964784L,
        -1.00000000000004092582630415831547637L,
        -0.999999999999934600951941089847743544L,
        -1.00000000000000654396874989191937317L,
        -0.999999999999999698757512863321320505L,
        -0.999999999999999996011105293688944503L,
        -1.00000000000000000185821504337850325L,
        -0.99999999999999999983971143614707759L,
        -1.00000000000000000000711558273921601L,
        -0.999999999999999999999992069688084904L,
        -0.999999999999999999999973441357813012L,
        -1.00000000000000000000000239515544094L,
        -0.999999999999999999999999881967211608L,
        -1.00000000000000000000000000229558469L,
        -1.00000000000000000000000000017456982L,
        -0.999999999999999999999999999978668328L,
        -1.00000000000000000000000000000124496L,
        -0.999999999999999999999999999999958088L,
        -1.00000000000000000000000000000000007L,
        -1.0000000000000000000000000000000001L,
        -0.999999999999999999999999999999999992L,
        -1.0L,
        -1.0L,
        -1.0L,
        -1.0L,
        -1.0L,
        -1.0L,
        -1.0L,
        -1.0L,
        -1.0L
};
        cplx<T> r{}, pw(T(1));
        for (const long double ck : c) {
            const cplx<T> t = static_cast<T>(ck) * pw;
            r += t;
            if (std::abs(t) <= eps<T>() / T(8) * std::abs(r)) break;
            pw *= s;
        }
        return r;
    }
    if (s.real() < T(0.5)) {
        // zeta(s) = 2^s pi^(s-1) sin(pi s/2) Gamma(1-s) zeta(1-s), in log form to avoid overflow.
        // The exponent is formed in double-T and sin(pi s/2) in split form, so neither
        // the phase |Im s| log|s| nor e^{pi |Im s|/2} costs accuracy or overflows.
        const cplx<T> one_minus_s = cplx<T>(T(1)) - s;
        const dd<T> ln2 = dd_const<T>(ddc::ln2), lnpi = dd_const<T>(ddc::log_pi);
        const dd<T> sm1 = two_sum(s.real(), T(-1));
        const cdd<T> lg1 = lgamma_dd(two_sum(T(1), -s.real()), -s.imag());
        const sinpi_exp<T> sp = sinpi_split(s / T(2));
        const cdd<T> e{ln2 * s.real() + lnpi * sm1 + lg1.re + sp.E,
                       ln2 * s.imag() + lnpi * s.imag() + lg1.im};
        return exp_split(e) * sp.P * riemann_zeta(one_minus_s);
    }
    const T t = std::abs(s.imag());
    if (t > T(1)) return zeta_em(s);
    const T target = T(std::numeric_limits<T>::digits) * std::log(T(2)) + T(3);
    const T rate = std::log(T(3) + detail::lm::sqrt(T(8)));
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
            // (k+1)^{-s} with log(k+1) in double-T: the phase Im(s) log(k+1) is large.
            const dd<T> lk = log_dd(dd<T>{T(k + 1), T()});
            const cplx<T> term = coeff * exp_split(cdd<T>{-(lk * s.real()), -(lk * s.imag())});
            eta += (k % 2 == 0) ? term : -term;
        }
    }
    // 1 - 2^(1-s) = -expm1((1-s) ln 2), accurate near s = 1.
    cplx<T> denom;
    if (std::abs(s - T(1)) < T(1)) {
        denom = -expm1((cplx<T>(T(1)) - s) * std::log(T(2)));
    } else {
        const dd<T> ln2 = dd_const<T>(ddc::ln2);
        denom = T(1) - exp_split(cdd<T>{ln2 * two_sum(T(1), -s.real()).hi + ln2 * two_sum(T(1), -s.real()).lo,
                                        -(ln2 * s.imag())});
    }
    return cdiv(eta, denom);
}

// ---- Ascending series for J and I, evaluated in double-T -----------------------------
// J_nu(z) = (z/2)^nu / Gamma(nu+1) * S,  S = sum_k (-+ z^2/4)^k / (k! (nu+1)_k)   (DLMF 10.2.2, 10.25.2)
// Used for |z| <= 16.  S is summed in double-T, so the cancellation between its terms (at
// most about e^|z|) costs nothing at T precision: the result has a *relative*
// error of a few ulp, including next to the zeros of J.  The prefactor is formed as
// exp(nu log(z/2) - log Gamma(nu+1)) with the exponent in double-T; for nu + 1 < 0 the
// imaginary part of the loggamma branch supplies the sign of Gamma.  nu must not be a
// negative integer (callers use J_{-n} = (-1)^n J_n, I_{-n} = I_n) and z != 0.
template <class T>
[[nodiscard]] inline bool bessel_series_ok(T nu, T az) noexcept {
    (void)nu;
    return az <= T(16);
}
template <class T>
[[nodiscard]] inline cplx<T> bessel_series(T nu, const cplx<T>& z, bool modified) noexcept {
    const dd<T> x{z.real() / T(2), T()}, y{z.imag() / T(2), T()};
    cdd<T> h2{two_prod(x.hi, x.hi) - two_prod(y.hi, y.hi), two_prod(x.hi, y.hi) * T(2)};
    if (!modified) h2 = {-h2.re, -h2.im};
    const T peak = std::abs(z) * std::abs(z) / T(4);
    const T tol = eps<T>() * eps<T>() / T(8);
    cdd<T> term{{T(1), T()}, {}}, sum = term;
    T maxterm = T(1);
    for (int k = 1; k < 100000; ++k) {
        const dd<T> inv = dd<T>{T(1), T()} / (two_sum(nu, T(k)) * T(k));
        term = term * h2;
        term = {term.re * inv, term.im * inv};
        sum = sum + term;
        const T at = std::max(std::abs(term.re.hi), std::abs(term.im.hi));
        maxterm = std::max(maxterm, at);
        if (T(k) * (T(k) + nu) > peak &&
            at <= tol * std::max(std::max(std::abs(sum.re.hi), std::abs(sum.im.hi)), eps<T>() * maxterm))
            break;
    }
    if (nu == T()) return to_cplx(sum);
    const T lmax = std::log(std::numeric_limits<T>::max()) / T(2);
    const T g = std::abs(nu) <= T(150) ? detail::lm::tgamma(nu + T(1)) : T();
    if (std::isfinite(g) && g != T() && std::abs(nu * std::log(std::abs(z) / T(2))) < lmax) {
        // (z/2)^nu = (z/2)^n (z/2)^f, n = floor(nu): the integer power by binary powering in
        // double-T, the fractional one (|f log(z/2)| small) in T; 1/Gamma(nu+1) from the C
        // library (correctly signed for nu + 1 < 0).
        const T fl = std::floor(nu), f = nu - fl;
        long n = static_cast<long>(std::abs(fl));
        cdd<T> b{x, y}, pw{{T(1), T()}, {}};
        while (n) {
            if (n & 1) pw = pw * b;
            b = b * b;
            n >>= 1;
        }
        cplx<T> pre = to_cplx(pw);
        if (fl < T()) pre = cdiv(T(1), pre);
        if (f != T()) pre *= detail::lm::exp(f * std::log(cplx<T>(x.hi, y.hi)));
        return (pre / g) * to_cplx(sum);
    }
    // Large order: exp(nu log(z/2) - log Gamma(nu+1)) with the exponent in double-T.
    const cdd<T> lz = clog_dd(x, y);
    const cdd<T> lg = lgamma_dd(two_sum(nu, T(1)), T());
    return exp_split(cdd<T>{lz.re * nu - lg.re, lz.im * nu - lg.im}) * to_cplx(sum);
}

template <class T>
struct bessel_all;
template <class T>
[[nodiscard]] inline bessel_all<T> bessel_jyh(T nu, const cplx<T>& z);
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
    if (bessel_series_ok(T(n), az)) return bessel_series(T(n), z, false);
    if (az >= T(16) && T(n) * T(n) < az) return bessel_jyh(T(n), z).j;  // Hankel asymptotics
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
    return cdiv(fn * detail::lm::cos(z), s2);
}

// std::acos with an exact pi/2 at zero (libc++ 19 returns a double-precision pi/2
// for complex<long double> there).
// asin and acos by the algorithm of Hull, Fairgrieve and Tang (ACM TOMS 23, 1997),
// accurate to a few ulp everywhere, including near the branch points +-1 and for huge
// or tiny arguments.  (libc++'s complex<long double> versions are only good to about
// double precision.)  The kernel works on x, y >= 0; the signs are restored with the
// C99 symmetries, so signed zeros select the side of the cuts.
template <class T>
inline void asin_acos_kernel(T x, T y, T& asin_re, T& acos_re, T& im) noexcept {
    const T safe_max = detail::lm::sqrt(std::numeric_limits<T>::max()) / T(8);
    const T safe_min = detail::lm::sqrt(std::numeric_limits<T>::min()) * T(4);
    if (y == T() && x <= T(1)) {
        asin_re = detail::lm::asin(x);
        acos_re = std::acos(x);
        im = T();
        return;
    }
    if (x > safe_max || y > safe_max) {
        // asin z ~ -i log(2 i z): error O(1/|z|^2), far below eps here.
        asin_re = std::atan2(x, y);
        acos_re = std::atan2(y, x);
        const T m = std::max(x, y), n = std::min(x, y), q = n / m;
        im = std::log(m) + std::log(T(2)) + T(0.5) * std::log1p(q * q);
        return;
    }
    if (x < safe_min && y < safe_min) {
        asin_re = x;
        acos_re = pi<T> / T(2) - x;
        im = y;
        return;
    }
    const T a_cross = T(1.5), b_cross = T(0.6417);
    const T xp1 = x + T(1), xm1 = x - T(1);
    const T r = std::hypot(xp1, y), s = std::hypot(xm1, y);
    const T A = (r + s) / T(2), B = x / A;
    // Real parts.
    if (B <= b_cross) {
        asin_re = detail::lm::asin(B);
        acos_re = std::acos(B);
    } else if (x <= T(1)) {
        const T d = detail::lm::sqrt(T(0.5) * (A + x) * (y * (y / (r + xp1)) + (s - xm1)));
        asin_re = std::atan(x / d);
        acos_re = std::atan(d / x);
    } else {
        const T h = T(0.5) * ((A + x) / (r + xp1) + (A + x) / (s + xm1));
        const T d = y * detail::lm::sqrt(h);
        asin_re = std::atan(x / d);
        acos_re = std::atan(d / x);
    }
    // Imaginary part (of asin; acos has the opposite sign).
    if (A <= a_cross) {
        T am1, sq;
        if (x < T(1)) {
            // sqrt(Am1) = y sqrt((1/(r+xp1) + 1/(s-xm1))/2), formed without y^2.
            const T sa = y * detail::lm::sqrt(T(0.5) * (T(1) / (r + xp1) + T(1) / (s - xm1)));
            am1 = sa * sa;
            sq = sa * detail::lm::sqrt(A + T(1));
        } else {
            am1 = T(0.5) * (y * (y / (r + xp1)) + (s + xm1));
            sq = detail::lm::sqrt(am1 * (A + T(1)));
        }
        im = std::log1p(am1 + sq);
    } else {
        im = std::log(A + detail::lm::sqrt(A * A - T(1)));
    }
}
template <class T>
[[nodiscard]] inline cplx<T> asin(const cplx<T>& z) noexcept {
    if (!std::isfinite(z.real()) || !std::isfinite(z.imag())) return std::asin(z);
    T ar, cr, im;
    asin_acos_kernel(std::abs(z.real()), std::abs(z.imag()), ar, cr, im);
    return {std::copysign(ar, z.real()), std::copysign(im, z.imag())};
}
template <class T>
[[nodiscard]] inline cplx<T> acos(const cplx<T>& z) noexcept {
    if (!std::isfinite(z.real()) || !std::isfinite(z.imag())) return std::acos(z);
    T ar, cr, im;
    asin_acos_kernel(std::abs(z.real()), std::abs(z.imag()), ar, cr, im);
    return {std::signbit(z.real()) ? pi<T> - cr : cr, std::signbit(z.imag()) ? im : -im};
}
// asinh z = -i asin(i z);  acosh z = +-i acos z with Re >= 0.
template <class T>
[[nodiscard]] inline cplx<T> asinh(const cplx<T>& z) noexcept {
    const cplx<T> w = sf::asin(cplx<T>(-z.imag(), z.real()));
    return {w.imag(), -w.real()};
}
template <class T>
[[nodiscard]] inline cplx<T> acosh(const cplx<T>& z) noexcept {
    if (!std::isfinite(z.real()) || !std::isfinite(z.imag())) return std::acosh(z);
    const cplx<T> w = sf::acos(z);
    return {std::abs(w.imag()), std::copysign(w.real(), z.imag())};
}

// atanh and atan with the C99/std::complex branch cuts, written in terms of log1p
// and cdiv (libc++'s complex<long double> versions lose precision through division).
//   atanh z = log1p(4x / ((1 - x)^2 + y^2)) / 4                 (real part, x >= 0)
//           + i atan2(2y, (1 - x)(1 + x) - y^2) / 2             for |x|, |y| < 2
//   atanh z = atanh(1/z) + i copysign(pi/2, Im z)              otherwise
//   atan  z = -i atanh(i z)
// The real part of (log1p z - log1p(-z)) / 2 cancelled for small Re z (4.5e15 eps in the
// vector part of atan(x + tiny v)), and through 1/z next to z = +-1 both parts inherited the
// rounding of the division in 1 - Re(1/z) (8.5e4 eps at z = 1 - 2.5e-14 + 3.2e-7 i; 1.4e4 eps
// at 1 + 7.2e-11 + 5.3e-10 i) before the FMA audit.  Here 1 - x is exact next to x = 1, and
// the products are relatively accurate whether or not they are contracted into fmas.
template <class T>
[[nodiscard]] inline cplx<T> atanh(const cplx<T>& z) noexcept {
    if (std::isnan(z.real()) || std::isnan(z.imag())) return std::atanh(z);
    const T x = z.real(), y = z.imag();
    if (std::abs(x) < T(2) && std::abs(y) < T(2)) {
        if (std::signbit(x)) return -sf::atanh(-z);  // atanh is odd
        const T xm = T(1) - x, d = xm * xm + y * y;
        return {std::log1p(T(4) * x / d) / T(4), std::atan2(T(2) * y, xm * (T(1) + x) - y * y) / T(2)};
    }
    const cplx<T> w = sf::atanh(cdiv(T(1), z));
    return {w.real(), w.imag() + std::copysign(pi<T> / T(2), z.imag())};
}
template <class T>
[[nodiscard]] inline cplx<T> atan(const cplx<T>& z) noexcept {
    const cplx<T> w = sf::atanh(cplx<T>(-z.imag(), z.real()));
    return {w.imag(), -w.real()};
}

// =====================================================================================
// Extended special functions: complex kernels
// =====================================================================================

template <class T>
[[nodiscard]] inline cplx<T> mul_i(const cplx<T>& z) noexcept { return {-z.imag(), z.real()}; }
template <class T>
[[nodiscard]] inline cplx<T> mul_mi(const cplx<T>& z) noexcept { return {z.imag(), -z.real()}; }
// exp(i pi x) for real x, exact at integers and half-integers.
template <class T>
[[nodiscard]] inline cplx<T> expi_pi(T x) noexcept { return {cospi(x), sinpi(x)}; }
template <class T>
[[nodiscard]] inline bool on_lower_side(const cplx<T>& z) noexcept {
    return z.imag() < T() || (z.imag() == T() && std::signbit(z.imag()));
}
template <class T>
[[nodiscard]] inline cplx<T> ipow(cplx<T> z, unsigned n) noexcept {
    cplx<T> r(T(1));
    for (; n; n >>= 1) {
        if (n & 1u) r *= z;
        z *= z;
    }
    return r;
}
template <class T>
[[nodiscard]] inline bool is_finite(const cplx<T>& z) noexcept {
    return std::isfinite(z.real()) && std::isfinite(z.imag());
}

// ---- double-exponential quadrature on (0, inf) ---------------------------------------
// Ooura-Mori transformation u = exp(t - exp(-t)) for integrands that decay
// exponentially and may have an integrable algebraic singularity at 0.  The step is
// halved until two levels agree to 32 eps relative to the sum of |terms|.
template <class T, class F>
[[nodiscard]] inline cplx<T> de_half_infinite(F&& f) {
    const T e = eps<T>();
    auto term = [&](T t) -> cplx<T> {
        const T et = detail::lm::exp(-t);
        const T u = detail::lm::exp(t - et);
        if (!(u > T()) || !std::isfinite(u)) return {};
        const cplx<T> v = f(u);
        return v * ((T(1) + et) * u);
    };
    T h = T(0.25);
    const T lo0 = T(-7), hi0 = T(12);
    const int n0 = static_cast<int>((hi0 - lo0) / h);
    std::vector<cplx<T>> first(static_cast<std::size_t>(n0) + 1);
    T maxmag{};
    for (int k = 0; k <= n0; ++k) {
        first[static_cast<std::size_t>(k)] = term(lo0 + T(k) * h);
        maxmag = std::max(maxmag, std::abs(first[static_cast<std::size_t>(k)]));
    }
    if (!std::isfinite(maxmag)) return {nan<T>(), nan<T>()};
    if (maxmag == T()) return {};
    int kl = 0, kh = n0;
    const T cut = maxmag * e * T(1e-6);
    while (kl < n0 && std::abs(first[static_cast<std::size_t>(kl)]) < cut) ++kl;
    while (kh > 0 && std::abs(first[static_cast<std::size_t>(kh)]) < cut) --kh;
    kl = std::max(0, kl - 1);
    kh = std::min(n0, kh + 1);
    const T lo = lo0 + T(kl) * h, hi = lo0 + T(kh) * h;
    cplx<T> sum{};
    T abssum{};
    for (int k = kl; k <= kh; ++k) {
        sum += first[static_cast<std::size_t>(k)];
        abssum += std::abs(first[static_cast<std::size_t>(k)]);
    }
    cplx<T> S = sum * h;
    for (int level = 1; level <= 11; ++level) {
        h /= T(2);
        cplx<T> add{};
        for (T t = lo + h; t < hi; t += T(2) * h) {
            const cplx<T> v = term(t);
            add += v;
            abssum += std::abs(v);
        }
        sum += add;
        const cplx<T> Snew = sum * h;
        const T scale = abssum * h;
        const bool done = std::abs(Snew - S) <= T(32) * e * scale;
        S = Snew;
        if (done && level >= 2) break;
    }
    return S;
}

// ---- Bessel functions of real order -------------------------------------------------
// K_mu(w), mu >= 0, |ph w| <= pi/2, w != 0:
//   K_mu(w) = sqrt(pi/(2w)) e^{-w} / Gamma(mu+1/2) * int_0^inf e^{-u} u^{mu-1/2} (1 + u/(2w))^{mu-1/2} du
// (DLMF 10.32.8 with the contour rotated so that the exponential is real).  The integrand
// does not oscillate, so the double-exponential rule converges for every |w|.
template <class T>
[[nodiscard]] inline cplx<T> bessel_k_laplace(T mu, const cplx<T>& w) {
    const T a = mu - T(0.5);
    const T lg = detail::lm::lgamma(mu + T(0.5));
    const cplx<T> inv2w = cdiv(T(0.5), w);
    const cplx<T> I = de_half_infinite<T>([&](T u) -> cplx<T> {
        const cplx<T> ex = cplx<T>(-u - lg) + a * (cplx<T>(std::log(u)) + log1p(u * inv2w));
        return detail::lm::exp(ex);
    });
    return detail::lm::sqrt(cdiv(pi<T>, T(2) * w)) * detail::lm::exp(-w) * I;
}

// Hankel's asymptotic expansion (DLMF 10.40.2)
//   K_mu(w) ~ sqrt(pi/(2w)) e^{-w} sum_k a_k(mu) / w^k,  a_k = prod_{j<=k} (4mu^2 - (2j-1)^2) / (8j),
// valid for |ph w| < 3pi/2 (the sign of a zero imaginary part of w selects ph w = +-pi).
// The error is below the first omitted term.  Returns false unless the terms decrease
// monotonically to below eps/4 of the sum (so it is used only where it is accurate);
// for half-integer mu the series terminates and is exact.
template <class T>
[[nodiscard]] inline bool bessel_k_asym(T mu, const cplx<T>& w, cplx<T>& k) {
    const T m4 = T(4) * mu * mu;
    const cplx<T> inv = cdiv(T(1), T(8) * w);
    cplx<T> term(T(1)), sum(T(1));
    T prev = T(1);
    for (int j = 1; j < 400; ++j) {
        const T odd = T(2 * j - 1);
        term *= ((m4 - odd * odd) / T(j)) * inv;
        const T at = std::abs(term);
        sum += term;
        if (at <= eps<T>() / T(4) * std::abs(sum)) {
            // sqrt(pi/2)/sqrt(w), not sqrt(pi/(2w)): the division would flip the sign of a
            // zero imaginary part and lose the ph w = -pi side used by I_mu.
            k = cdiv(detail::lm::sqrt(pi<T> / T(2)) * detail::lm::exp(-w), detail::lm::sqrt(w)) * sum;
            return true;
        }
        if (!(at < prev)) return false;
        prev = at;
    }
    return false;
}
// K_mu(w) for Re w >= 0: the asymptotic expansion for large |w|, otherwise the
// double-exponential Laplace integral.
template <class T>
[[nodiscard]] inline cplx<T> bessel_k_one(T mu, const cplx<T>& w) {
    cplx<T> k;
    if (std::abs(w) >= T(8) && bessel_k_asym(mu, w, k)) return k;
    return bessel_k_laplace(mu, w);
}
// K_mu and K_{mu+1}, Re w >= 0.  For mu >= 2, K is computed at the orders mu0 and mu0+1
// (mu0 = frac(mu)) and recurred upward, K_{n+1} = K_{n-1} + (2n/w) K_n; K is the dominant
// solution for increasing order on Re w >= 0, so the recurrence is stable (Temme).
template <class T>
struct kpair { cplx<T> k0, k1; };
// K_a and K_{a+1} for |a| <= 1/2 and small |w| (|w| <= 2, Re w >= 0) by Temme's series
// (Temme 1975; Numerical Recipes 6.7), with 1/Gamma(1 +- a) and the auxiliary
// gamma1(a), gamma2(a) taken from the Taylor series of 1/Gamma(1+x).
template <class T>
[[nodiscard]] inline kpair<T> bessel_k_temme(T a, const cplx<T>& w) noexcept {
    static constexpr long double c[] = {
    1.0L, 0.577215664901532860606512090082402431L, -0.655878071520253881077019515145390481L,
    -0.0420026350340952355290039348754298187L, 0.166538611382291489501700795102105236L,
    -0.0421977345555443367482083012891873913L, -0.00962197152787697356211492167234819898L,
    0.00721894324666309954239501034044657271L, -0.00116516759185906511211397108401838867L,
    -0.000215241674114950972815729963053647806L, 0.000128050282388116186153198626328164323L,
    -0.0000201348547807882386556893914210218184L, -0.00000125049348214267065734535947383309224L,
    0.00000113302723198169588237412962033074494L, -0.000000205633841697760710345015413002057284L,
    0.00000000611609510448141581786249868285534287L, 0.0000000050020076444692229300556650480599913L,
    -0.00000000118127457048702014458812656543650558L, 0.000000000104342671169110051049154033231225019L,
    7.78226343990507125404993731136077723e-12L, -3.69680561864220570818781587808576624e-12L,
    5.1003702874544759790154813228632318e-13L, -2.05832605356650678322242954485523742e-14L,
    -5.34812253942301798237001731872793995e-15L, 1.22677862823826079015889384662242243e-15L,
    -1.18125930169745876951376458684229783e-16L, 1.18669225475160033257977724292867407e-18L,
    1.41238065531803178155580394756670904e-18L, -2.29874568443537020659247858063369926e-19L,
    1.71440632192733743338396337026725707e-20L, 1.33735173049369311486478139512226802e-22L,
    -2.0542335517666727893250253513557338e-22L, 2.73603004860799984483150990433098201e-23L};
    T gp{}, gm{}, g1{}, g2{}, ak = T(1);
    for (int k = 0; k < static_cast<int>(sizeof(c) / sizeof(c[0])); ++k) {
        const T ck = static_cast<T>(c[k]);
        gp += ck * ak;
        gm += (k % 2 ? -ck : ck) * ak;
        if (k % 2) g1 -= ck * ak / (a == T() ? T(1) : a);
        else g2 += ck * ak;
        ak *= a;
    }
    if (a == T()) g1 = -static_cast<T>(c[1]);
    const cplx<T> x2 = w / T(2);
    const cplx<T> d = -std::log(x2);
    const cplx<T> e = a * d;
    const T fact = a == T() ? T(1) : a * pi<T> / sinpi(a);
    const cplx<T> fact2 = std::abs(e) < eps<T>() ? cplx<T>(T(1)) : cdiv(detail::lm::sinh(e), e);
    cplx<T> ff = fact * (g1 * detail::lm::cosh(e) + g2 * fact2 * d);
    cplx<T> sum = ff;
    const cplx<T> ee = detail::lm::exp(e);
    cplx<T> p = T(0.5) * ee / gp, q = cdiv(T(0.5), ee * gm);  // cdiv: libc++ divides long double complex in double
    cplx<T> cc(T(1));
    const cplx<T> x22 = x2 * x2;
    cplx<T> sum1 = p, c0{}, c1{};
    for (int i = 1; i < 500; ++i) {
        const T ti = T(i);
        ff = (ti * ff + p + q) / (ti * ti - a * a);
        cc *= x22 / ti;
        p /= (ti - a);
        q /= (ti + a);
        const cplx<T> del = cc * ff, del1 = cc * (p - ti * ff);
        const cplx<T> y0 = del - c0, t0 = sum + y0;  // compensated sums
        c0 = (t0 - sum) - y0;
        sum = t0;
        const cplx<T> y1 = del1 - c1, t1 = sum1 + y1;
        c1 = (t1 - sum1) - y1;
        sum1 = t1;
        if (std::abs(del) < eps<T>() / T(4) * std::abs(sum) && std::abs(del1) < eps<T>() / T(4) * std::abs(sum1)) break;
    }
    return {sum, sum1 * cdiv(T(2), w)};
}
// K_a and K_{a+1}, |a| <= 1/2, Re w >= 0, by Steed's method on Temme's continued fraction
// CF2 (Numerical Recipes 6.7, Thompson & Barnett 1987); used for 2 < |w| < 8.  Sets ok to
// false if it has not converged.
template <class T>
[[nodiscard]] inline kpair<T> bessel_k_cf2(T a, const cplx<T>& w, bool& ok) noexcept {
    const T a1 = T(0.25) - a * a;
    cplx<T> b = T(2) * (T(1) + w), d = cdiv(T(1), b), h = d, delh = d;
    cplx<T> q1{}, q2(T(1)), q(a1), c(a1), s = T(1) + q * delh, hc{}, sc{};
    T aa = -a1;
    ok = false;
    for (int i = 1; i < 20000; ++i) {
        aa -= T(2 * i);
        c = -aa * c / T(i + 1);
        const cplx<T> qnew = (q1 - b * q2) / aa;
        q1 = q2;
        q2 = qnew;
        q += c * qnew;
        b += T(2);
        d = cdiv(T(1), b + aa * d);
        delh = (b * d - T(1)) * delh;
        // compensated (Kahan) sums: the terms alternate in phase and many are needed near
        // the imaginary axis
        const cplx<T> yh = delh - hc, th = h + yh;
        hc = (th - h) - yh;
        h = th;
        const cplx<T> dels = q * delh;
        const cplx<T> ys = dels - sc, ts = s + ys;
        sc = (ts - s) - ys;
        s = ts;
        if (std::abs(dels) < eps<T>() / T(4) * std::abs(s) && std::abs(delh) < eps<T>() / T(4) * std::abs(h)) {
            ok = true;
            break;
        }
    }
    h *= a1;
    const cplx<T> k0 = cdiv(detail::lm::sqrt(pi<T> / T(2)) * detail::lm::exp(-w), detail::lm::sqrt(w) * s);
    const cplx<T> k1 = k0 * cdiv(a + w + T(0.5) - h, w);
    return {k0, k1};
}
template <class T>
[[nodiscard]] inline kpair<T> bessel_k_pair(T mu, const cplx<T>& w) {
    const T aw = std::abs(w);
    bool ok = true;
    kpair<T> r{};
    const T nr = std::nearbyint(mu), ar = mu - nr;
    // Hankel's expansion where it converges (|w| >= 8), otherwise the series / CF2 at the base
    // order ar and upward recurrence; the Laplace integral only if CF2 fails to converge.
    bool small = aw < T(8);
    if (!small) small = !(bessel_k_asym(ar, w, r.k0) && bessel_k_asym(ar + T(1), w, r.k1));
    if (small) {
        if (std::abs(ar) == T(0.5)) {
            // half-integer order: K_{1/2}(w) = sqrt(pi/(2w)) e^{-w}, K_{3/2} = K_{1/2} (1 + 1/w)
            const cplx<T> k = cdiv(detail::lm::sqrt(pi<T> / T(2)) * detail::lm::exp(-w), detail::lm::sqrt(w));
            r = ar > T() ? kpair<T>{k, k + cdiv(k, w)} : kpair<T>{k, k};  // K_{-1/2} = K_{1/2}
        } else if (aw <= T(2)) r = bessel_k_temme(ar, w);
        else r = bessel_k_cf2(ar, w, ok);
    }
    if (ok) {
        const cplx<T> tw = cdiv(T(2), w);
        for (long j = 1; j <= static_cast<long>(nr); ++j) {
            const cplx<T> c = r.k0 + (ar + T(j)) * tw * r.k1;
            r = {r.k1, c};
        }
        return r;
    }
    if (mu < T(2)) return {bessel_k_one(mu, w), bessel_k_one(mu + T(1), w)};
    const T fl = std::floor(mu), m0 = mu - fl;
    const long n = static_cast<long>(fl);
    cplx<T> a = bessel_k_one(m0, w), b = bessel_k_one(m0 + T(1), w);
    const cplx<T> two_over_w = cdiv(T(2), w);
    for (long j = 1; j <= n; ++j) {
        const cplx<T> c = a + (m0 + T(j)) * two_over_w * b;
        a = b;
        b = c;
    }
    return {a, b};
}

// Ratio Z_{mu+1}/Z_mu of the minimal solution by the continued fraction
//   1/(b1 + s/(b2 + s/(b3 + ...))), b_k = 2(mu+k)/z,
// with s = -1 for J and s = +1 for I (modified Lentz).
template <class T>
[[nodiscard]] inline cplx<T> bessel_ratio(T mu, const cplx<T>& z, T s) {
    const T tiny = std::numeric_limits<T>::min() * T(1e10);
    const cplx<T> two_over_z = cdiv(T(2), z);
    cplx<T> f = (mu + T(1)) * two_over_z;
    if (std::abs(f) < tiny) f = tiny;
    cplx<T> C = f, D{};
    const long kmax = 1000000 + static_cast<long>(4 * std::abs(z));
    for (long k = 2; k < kmax; ++k) {
        const cplx<T> b = (mu + T(k)) * two_over_z;
        D = b + s * D;
        if (std::abs(D) < tiny) D = tiny;
        D = cdiv(T(1), D);
        C = b + cdiv(s, C);
        if (std::abs(C) < tiny) C = tiny;
        const cplx<T> delta = C * D;
        f *= delta;
        if (std::abs(delta - T(1)) <= eps<T>()) break;
    }
    return cdiv(T(1), f);
}

template <class T>
struct bessel_jh { cplx<T> j, h1; };

// J_mu(z) and H1_mu(z) for z in the closed first quadrant, z != 0, mu >= 0:
//   H1_mu(z) = (2/(pi i)) e^{-i mu pi/2} K_mu(-iz)                 (DLMF 10.27.8)
//   J_mu(z)  = 2i / (pi z (rho H1_mu - H1_{mu+1})),  rho = J_{mu+1}/J_mu   (Wronskian)
template <class T>
[[nodiscard]] inline bessel_jh<T> bessel_core_q1(T mu, const cplx<T>& z) {
    const cplx<T> w = mul_mi(z);
    const cplx<T> c(T(), -T(2) / pi<T>);
    const kpair<T> kp = bessel_k_pair(mu, w);
    const cplx<T> h0 = c * expi_pi(-mu / T(2)) * kp.k0;
    const cplx<T> h1 = c * expi_pi(-(mu + T(1)) / T(2)) * kp.k1;
    // Large |z|: J = (H1 + H2)/2 with H2_mu(z) = (2i/pi) e^{i mu pi/2} K_mu(iz) from the
    // asymptotic expansion (ph(iz) in [pi/2, pi]); this avoids the O(|z|) continued fraction.
    cplx<T> k2;
    if (std::abs(z) >= T(16) && bessel_k_asym(mu, mul_i(z), k2)) {
        const cplx<T> h2 = cplx<T>(T(), T(2) / pi<T>) * expi_pi(mu / T(2)) * k2;
        return {(h0 + h2) / T(2), h0};
    }
    // Near overflow (large mu, tiny z) H1_{mu+1} can overflow while H1_mu does not; then
    // |J_mu| ~ 2/(pi |z H1_{mu+1}|) is below the smallest normal and is returned as 0
    // (instead of NaN from inf - inf).
    if (is_finite(h0) && !is_finite(h1)) return {cplx<T>(), h0};
    const cplx<T> rho = bessel_ratio(mu, z, T(-1));
    const cplx<T> j = cdiv(cplx<T>(T(), T(2) / pi<T>), z * (rho * h0 - h1));
    return {j, h0};
}

template <class T>
struct bessel_all { cplx<T> j, y, h1, h2; };

// J, Y, H1, H2 of real order nu at any z != 0 (principal branches, cut along the
// negative real axis; the side is taken from the sign of Im z, including signed zero).
template <class T>
[[nodiscard]] inline bessel_all<T> bessel_jyh(T nu, const cplx<T>& z) {
    const T mu = std::abs(nu);
    const bool lower = on_lower_side(z);
    const cplx<T> zz = lower ? std::conj(z) : z;
    cplx<T> j, h1;
    if (zz.real() >= T()) {
        const auto c = bessel_core_q1(mu, zz);
        j = c.j;
        h1 = c.h1;
    } else {
        // z = conj(w) e^{i pi}: J(z) = e^{i mu pi} conj J(w), H1(z) = -e^{-i mu pi} conj H1(w).
        const auto c = bessel_core_q1(mu, cplx<T>(-zz.real(), zz.imag()));
        j = expi_pi(mu) * std::conj(c.j);
        h1 = -expi_pi(-mu) * std::conj(c.h1);
    }
    bessel_all<T> r{j, mul_i(j - h1), h1, T(2) * j - h1};
    if (lower) r = {std::conj(r.j), std::conj(r.y), std::conj(r.h2), std::conj(r.h1)};
    if (nu < T()) {
        const T s = sinpi(mu), c = cospi(mu);
        r = {c * r.j - s * r.y, s * r.j + c * r.y, expi_pi(mu) * r.h1, expi_pi(-mu) * r.h2};
    }
    return r;
}

template <class T>
[[nodiscard]] inline bool is_integer(T x) noexcept { return std::isfinite(x) && x == std::floor(x); }

template <class T>
[[nodiscard]] inline cplx<T> bessel_j_nu(T nu, const cplx<T>& z) {
    if (!is_finite(z)) return {nan<T>(), nan<T>()};
    if (z == cplx<T>()) {
        if (nu == T()) return {T(1), T()};
        if (nu > T() || is_integer(nu)) return {};
        return {std::numeric_limits<T>::infinity(), T()};
    }
    if (is_integer(nu) && nu < T()) {
        const cplx<T> r = bessel_j_nu(-nu, z);
        return std::fmod(nu, T(2)) == T() ? r : -r;
    }
    if (bessel_series_ok(nu, std::abs(z))) return bessel_series(nu, z, false);
    return bessel_jyh(nu, z).j;
}
template <class T>
[[nodiscard]] inline cplx<T> bessel_y_nu(T nu, const cplx<T>& z) {
    if (!is_finite(z)) return {nan<T>(), nan<T>()};
    if (z == cplx<T>()) return {-std::numeric_limits<T>::infinity(), T()};
    return bessel_jyh(nu, z).y;
}
template <class T>
[[nodiscard]] inline cplx<T> hankel_1(T nu, const cplx<T>& z) {
    if (z == cplx<T>()) return {nan<T>(), -std::numeric_limits<T>::infinity()};
    return bessel_jyh(nu, z).h1;
}
template <class T>
[[nodiscard]] inline cplx<T> hankel_2(T nu, const cplx<T>& z) {
    if (z == cplx<T>()) return {nan<T>(), std::numeric_limits<T>::infinity()};
    return bessel_jyh(nu, z).h2;
}

// K_nu(w) for any w != 0 (cut along the negative real axis).
template <class T>
[[nodiscard]] inline cplx<T> bessel_k_nu(T nu, const cplx<T>& w) {
    const T mu = std::abs(nu);
    if (!is_finite(w)) return {nan<T>(), nan<T>()};
    if (w == cplx<T>()) return {std::numeric_limits<T>::infinity(), T()};
    if (on_lower_side(w)) return std::conj(bessel_k_nu(mu, std::conj(w)));
    if (w.real() >= T()) return bessel_k_pair(mu, w).k0;
    // K_mu(w) = -(pi i/2) e^{-i mu pi/2} H2_mu(-iw), -pi/2 < ph w <= pi   (DLMF 10.27.8)
    return cplx<T>(T(), -pi<T> / T(2)) * expi_pi(-mu / T(2)) * bessel_jyh(mu, mul_mi(w)).h2;
}

// I_mu(w), mu >= 0, Re w >= 0, w != 0, from the Wronskian I_mu K_{mu+1} + I_{mu+1} K_mu = 1/w.
template <class T>
[[nodiscard]] inline cplx<T> bessel_i_right(T mu, const cplx<T>& w) {
    // Large |w|: I = (K_mu(w e^{-i pi}) - e^{i mu pi} K_mu(w)) / (pi i)  (DLMF 10.34.2), both
    // K from the asymptotic expansion; Im w < 0 by conjugation.
    if (std::abs(w) >= T(16)) {
        const bool low = on_lower_side(w);
        const cplx<T> ww = low ? std::conj(w) : w;
        // -ww with ph in [-pi, -pi/2]: a zero imaginary part becomes -0
        const cplx<T> wm(-ww.real(), -ww.imag());
        cplx<T> ka, kb;
        if (bessel_k_asym(mu, wm, ka) && bessel_k_asym(mu, ww, kb)) {
            const cplx<T> r = mul_mi(ka - expi_pi(mu) * kb) / pi<T>;
            return low ? std::conj(r) : r;
        }
    }
    const cplx<T> r = bessel_ratio(mu, w, T(1));
    const kpair<T> kp = bessel_k_pair(mu, w);
    return cdiv(T(1), w * (kp.k1 + r * kp.k0));
}
template <class T>
[[nodiscard]] inline cplx<T> bessel_i_nu(T nu, const cplx<T>& w) {
    const T mu = std::abs(nu);
    if (!is_finite(w)) return {nan<T>(), nan<T>()};
    if (w == cplx<T>()) {
        if (nu == T()) return {T(1), T()};
        if (nu > T() || is_integer(nu)) return {};
        return {std::numeric_limits<T>::infinity(), T()};
    }
    if (bessel_series_ok(is_integer(nu) ? mu : nu, std::abs(w))) return bessel_series(is_integer(nu) ? mu : nu, w, true);
    cplx<T> i;
    if (w.real() >= T()) i = bessel_i_right(mu, w);
    else i = expi_pi(on_lower_side(w) ? -mu : mu) * bessel_i_right(mu, -w);  // I(w e^{+-i pi}) = e^{+-i mu pi} I(w)
    if (nu < T() && !is_integer(mu)) i += (T(2) / pi<T>) * sinpi(mu) * bessel_k_nu(mu, w);
    return i;
}

// Spherical Bessel functions: f_n(z) = sqrt(pi/(2z)) F_{n+1/2}(z).
template <class T>
[[nodiscard]] inline cplx<T> sph_factor(const cplx<T>& z) { return detail::lm::sqrt(cdiv(pi<T>, T(2) * z)); }
template <class T>
[[nodiscard]] inline cplx<T> sph_bessel_j(unsigned n, const cplx<T>& z) {
    if (z == cplx<T>()) return {n == 0 ? T(1) : T(), T()};
    if (bessel_series_ok(T(n) + T(0.5), std::abs(z))) return sph_factor(z) * bessel_series(T(n) + T(0.5), z, false);
    return sph_factor(z) * bessel_jyh(T(n) + T(0.5), z).j;
}
template <class T>
[[nodiscard]] inline cplx<T> sph_bessel_y(unsigned n, const cplx<T>& z) {
    if (z == cplx<T>()) return {-std::numeric_limits<T>::infinity(), T()};
    return sph_factor(z) * bessel_jyh(T(n) + T(0.5), z).y;
}
template <class T>
[[nodiscard]] inline cplx<T> sph_hankel(unsigned n, const cplx<T>& z, int kind) {
    if (z == cplx<T>()) return {nan<T>(), nan<T>()};
    const auto b = bessel_jyh(T(n) + T(0.5), z);
    return sph_factor(z) * (kind == 1 ? b.h1 : b.h2);
}

// ---- Airy functions ---------------------------------------------------------------------
// |ph z| <= 2pi/3: Ai = sqrt(z/3) K_{1/3}(zeta)/pi, Bi = sqrt(z/3) (2 I_{1/3}(zeta) + sqrt3/pi K_{1/3}(zeta)),
// otherwise with t = -z: Ai = sqrt(t)/2 (J_{1/3} - Y_{1/3}/sqrt3), Bi = -sqrt(t)/2 (J_{1/3}/sqrt3 + Y_{1/3}),
// zeta = (2/3) z^{3/2} (DLMF 9.6).
template <class T>
[[nodiscard]] inline cplx<T> airy(const cplx<T>& z, bool bi) {
    if (!is_finite(z)) return {nan<T>(), nan<T>()};
    if (z == cplx<T>()) {
        const T ai0 = static_cast<T>(0.355028053887817239260063186004183176L);
        const T bi0 = static_cast<T>(0.614926627446000735150922369093613554L);
        return {bi ? bi0 : ai0, T()};
    }
    if (on_lower_side(z) && z.imag() != T()) return std::conj(airy(std::conj(z), bi));
    const T third = T(1) / T(3), sqrt3 = detail::lm::sqrt(T(3));
    const T ph = std::arg(z);
    if (std::abs(ph) <= T(2) * pi<T> / T(3)) {
        const cplx<T> zeta = (T(2) / T(3)) * z * detail::lm::sqrt(z);
        const cplx<T> pre = detail::lm::sqrt(z / T(3));
        const cplx<T> k = bessel_k_nu(third, zeta);
        if (!bi) return pre * k / pi<T>;
        return pre * (T(2) * bessel_i_nu(third, zeta) + (sqrt3 / pi<T>) * k);
    }
    const cplx<T> t = -z;
    const cplx<T> zeta = (T(2) / T(3)) * t * detail::lm::sqrt(t);
    const auto b = bessel_jyh(third, zeta);
    const cplx<T> st = detail::lm::sqrt(t) / T(2);
    if (!bi) return st * (b.j - b.y / sqrt3);
    return -st * (b.j / sqrt3 + b.y);
}

// ---- Carlson symmetric elliptic integrals (complex arguments) ---------------------------
template <class T>
[[nodiscard]] inline cplx<T> carlson_rf(cplx<T> x, cplx<T> y, cplx<T> z) {
    const cplx<T> A0 = (x + y + z) / T(3);
    T Q = detail::lm::pow(T(3) * eps<T>(), T(-1) / T(6)) *
          std::max({std::abs(A0 - x), std::abs(A0 - y), std::abs(A0 - z)});
    cplx<T> A = A0;
    const cplx<T> x0 = x, y0 = y;
    T f = T(1);  // 4^{-m}
    for (int m = 0; m < 100 && !(f * Q < std::abs(A)); ++m) {
        const cplx<T> sx = detail::lm::sqrt(x), sy = detail::lm::sqrt(y), sz = detail::lm::sqrt(z);
        const cplx<T> lam = sx * sy + sx * sz + sy * sz;
        A = (A + lam) / T(4);
        x = (x + lam) / T(4);
        y = (y + lam) / T(4);
        z = (z + lam) / T(4);
        f /= T(4);
    }
    const cplx<T> X = cdiv((A0 - x0) * f, A), Y = cdiv((A0 - y0) * f, A), Z = -(X + Y);
    const cplx<T> E2 = X * Y - Z * Z, E3 = X * Y * Z;
    const cplx<T> s = T(1) - E2 / T(10) + E3 / T(14) + E2 * E2 / T(24) - T(3) * E2 * E3 / T(44);
    return cdiv(s, detail::lm::sqrt(A));
}
template <class T>
[[nodiscard]] inline cplx<T> carlson_rd(cplx<T> x, cplx<T> y, cplx<T> z) {
    const cplx<T> A0 = (x + y + T(3) * z) / T(5);
    T Q = detail::lm::pow(eps<T>() / T(4), T(-1) / T(6)) *
          std::max({std::abs(A0 - x), std::abs(A0 - y), std::abs(A0 - z)});
    cplx<T> A = A0, sum{};
    const cplx<T> x0 = x, y0 = y;
    T f = T(1);
    for (int m = 0; m < 100 && !(f * Q < std::abs(A)); ++m) {
        const cplx<T> sx = detail::lm::sqrt(x), sy = detail::lm::sqrt(y), sz = detail::lm::sqrt(z);
        const cplx<T> lam = sx * sy + sx * sz + sy * sz;
        sum += cdiv(cplx<T>(f), sz * (z + lam));
        A = (A + lam) / T(4);
        x = (x + lam) / T(4);
        y = (y + lam) / T(4);
        z = (z + lam) / T(4);
        f /= T(4);
    }
    const cplx<T> X = cdiv((A0 - x0) * f, A), Y = cdiv((A0 - y0) * f, A), Z = -(X + Y) / T(3);
    const cplx<T> E2 = X * Y - T(6) * Z * Z, E3 = (T(3) * X * Y - T(8) * Z * Z) * Z;
    const cplx<T> E4 = T(3) * (X * Y - Z * Z) * Z * Z, E5 = X * Y * Z * Z * Z;
    const cplx<T> s = T(1) - T(3) * E2 / T(14) + E3 / T(6) + T(9) * E2 * E2 / T(88) - T(3) * E4 / T(22) -
                      T(9) * E2 * E3 / T(52) + T(3) * E5 / T(26);
    return f * cdiv(s, A * detail::lm::sqrt(A)) + T(3) * sum;
}
// R_C(1, 1 + e)
template <class T>
[[nodiscard]] inline cplx<T> carlson_rc1(const cplx<T>& e) {
    return carlson_rf(cplx<T>(T(1)), T(1) + e, T(1) + e);
}
// R_J(x, y, z, p): x, y, z real non-negative (at most one zero), p complex in the cut plane
// (Carlson 1995, with the R_C(1, 1+e) form of the sum).
template <class T>
[[nodiscard]] inline cplx<T> carlson_rj(cplx<T> x, cplx<T> y, cplx<T> z, cplx<T> p) {
    const cplx<T> A0 = (x + y + z + T(2) * p) / T(5);
    const cplx<T> delta = (p - x) * (p - y) * (p - z);
    T Q = detail::lm::pow(eps<T>() / T(4), T(-1) / T(6)) *
          std::max({std::abs(A0 - x), std::abs(A0 - y), std::abs(A0 - z), std::abs(A0 - p)});
    cplx<T> A = A0, sum{};
    const cplx<T> x0 = x, y0 = y, z0 = z;
    T f = T(1);
    for (int m = 0; m < 100 && !(f * Q < std::abs(A)); ++m) {
        const cplx<T> sx = detail::lm::sqrt(x), sy = detail::lm::sqrt(y), sz = detail::lm::sqrt(z), sp = detail::lm::sqrt(p);
        const cplx<T> lam = sx * sy + sx * sz + sy * sz;
        const cplx<T> d = (sp + sx) * (sp + sy) * (sp + sz);
        const cplx<T> e = cdiv(delta * (f * f * f), d * d);
        sum += f * cdiv(carlson_rc1(e), d);
        A = (A + lam) / T(4);
        x = (x + lam) / T(4);
        y = (y + lam) / T(4);
        z = (z + lam) / T(4);
        p = (p + lam) / T(4);
        f /= T(4);
    }
    const cplx<T> X = cdiv((A0 - x0) * f, A), Y = cdiv((A0 - y0) * f, A), Z = cdiv((A0 - z0) * f, A);
    const cplx<T> P = -(X + Y + Z) / T(2);
    const cplx<T> E2 = X * Y + X * Z + Y * Z - T(3) * P * P;
    const cplx<T> E3 = X * Y * Z + T(2) * E2 * P + T(4) * P * P * P;
    const cplx<T> E4 = (T(2) * X * Y * Z + E2 * P + T(3) * P * P * P) * P;
    const cplx<T> E5 = X * Y * Z * P * P;
    const cplx<T> s = T(1) - T(3) * E2 / T(14) + E3 / T(6) + T(9) * E2 * E2 / T(88) - T(3) * E4 / T(22) -
                      T(9) * E2 * E3 / T(52) + T(3) * E5 / T(26);
    return f * cdiv(s, A * detail::lm::sqrt(A)) + T(6) * sum;
}

// p = 1 - nu c for real c > 0, keeping the side of the cut: a real nu > 1/c lies on the cut
// of R_J, and the upper side of nu maps to the lower side of p.
template <class T>
[[nodiscard]] inline cplx<T> rj_param(const cplx<T>& nu, T c) {
    cplx<T> p(T(1) - nu.real() * c, -nu.imag() * c);
    if (p.imag() == T() && p.real() < T()) p.imag(-std::numeric_limits<T>::min() * T(1e100));
    return p;
}
// The same with s = sin phi, c = cos phi for real phi: next to the pole nu sin^2 phi = 1,
// 1 - nu*s*s cancels, and how much it loses depends on whether it is contracted into an fma.
// Where |c| <= |s| (Re nu <= 2 near the pole) it is (1 - Re nu) + Re nu c^2, with 1 - Re nu
// exact and the small c^2 relatively accurate; otherwise s^2 and its product with Re nu are
// carried exactly (two_prod).
template <class T>
[[nodiscard]] inline cplx<T> rj_param_sin(const cplx<T>& nu, T s, T c) {
    const T nr = nu.real();
    T re;
    if (std::abs(c) <= std::abs(s)) {
        const dd<T> c2 = two_prod(c, c);
        re = (T(1) - nr) + (nr * c2.hi + nr * c2.lo);
    } else {
        const dd<T> s2 = two_prod(s, s), p1 = two_prod(nr, s2.hi);
        re = ((T(1) - p1.hi) - p1.lo) - nr * s2.lo;
    }
    cplx<T> p(re, -nu.imag() * (s * s));
    if (p.imag() == T() && p.real() < T()) p.imag(-std::numeric_limits<T>::min() * T(1e100));
    return p;
}
// k'^2 = 1 - k^2 as (1 - k)(1 + k): next to k = +-1 the factor 1 -+ k is exact (Sterbenz),
// while 1 - k*k loses the low bits of k*k (up to 7.3e5 eps in K and 7.6e5 eps in Pi next to
// k = +-1 before the FMA audit), and whether it does depends on FMA contraction
// (fma(-k, k, 1) is exact).
template <class T>
[[nodiscard]] inline cplx<T> one_minus_sq(const cplx<T>& k) {
    return (T(1) - k) * (T(1) + k);
}
template <class T>
[[nodiscard]] inline T one_minus_sq(T k) {
    return (T(1) - k) * (T(1) + k);
}
template <class T>
[[nodiscard]] inline cplx<T> comp_ellint_1k(const cplx<T>& k) {
    const cplx<T> y = one_minus_sq(k);
    if (y == cplx<T>()) return {std::numeric_limits<T>::infinity(), T()};
    return carlson_rf(cplx<T>(), y, cplx<T>(T(1)));
}
template <class T>
[[nodiscard]] inline cplx<T> comp_ellint_2k(const cplx<T>& k) {
    const cplx<T> y = one_minus_sq(k);
    if (y == cplx<T>()) return {T(1), T()};
    return carlson_rf(cplx<T>(), y, cplx<T>(T(1))) - k * k / T(3) * carlson_rd(cplx<T>(), y, cplx<T>(T(1)));
}
// Pi(nu, k) with std's convention int_0^{pi/2} dt / ((1 - nu sin^2 t) sqrt(1 - k^2 sin^2 t)).
template <class T>
[[nodiscard]] inline cplx<T> comp_ellint_3n(T k, const cplx<T>& nu) {
    const cplx<T> y(one_minus_sq(k));
    if (nu == cplx<T>(T(1))) return {std::numeric_limits<T>::infinity(), T()};
    return carlson_rf(cplx<T>(), y, cplx<T>(T(1))) +
           nu / T(3) * carlson_rj(cplx<T>(), y, cplx<T>(T(1)), rj_param(nu, T(1)));
}
// 1 - k^2 sin^2 phi as k'^2 + (k cos phi)^2: no cancellation next to k = 1, phi = pi/2, where
// 1 - k*k*s*s lost up to 5.4e11 eps in F (the amount depending on FMA contraction).
template <class T>
[[nodiscard]] inline cplx<T> ellint_delta2(T k, const cplx<T>& c) {
    const cplx<T> kc = k * c;
    return one_minus_sq(k) + kc * kc;
}
template <class T>
[[nodiscard]] inline T ellint_delta2(T k, T c) {
    const T kc = k * c;
    return one_minus_sq(k) + kc * kc;
}
// Reduce the amplitude: phi = phi0 + m pi with |Re phi0| <= pi/2.
template <class T>
[[nodiscard]] inline T amplitude_periods(const cplx<T>& phi) noexcept {
    return std::round(phi.real() / pi<T>);
}
// sin and cos of phi0 = phi - m pi as (-1)^m sin(phi), (-1)^m cos(phi): next to phi0 = +-pi/2
// the cosine is small, and cos(phi - m*pi) would inherit the rounding of m*pi and of the
// subtraction.
template <class T>
[[nodiscard]] inline T period_sign(T m) noexcept {
    return std::fmod(m, T(2)) == T() ? T(1) : T(-1);
}
template <class T>
[[nodiscard]] inline cplx<T> ellint_1p(T k, const cplx<T>& phi) {
    const T m = amplitude_periods(phi), sg = period_sign(m);
    const cplx<T> s = sg * detail::lm::sin(phi), c = sg * detail::lm::cos(phi);
    cplx<T> r = s * carlson_rf(c * c, ellint_delta2(k, c), cplx<T>(T(1)));
    if (m != T()) r += T(2) * m * comp_ellint_1k(cplx<T>(k));
    return r;
}
template <class T>
[[nodiscard]] inline cplx<T> ellint_2p(T k, const cplx<T>& phi) {
    const T m = amplitude_periods(phi), sg = period_sign(m);
    const cplx<T> s = sg * detail::lm::sin(phi), c = sg * detail::lm::cos(phi);
    const cplx<T> x = c * c, y = ellint_delta2(k, c);
    cplx<T> r = s * carlson_rf(x, y, cplx<T>(T(1))) -
                (k * k / T(3)) * s * s * s * carlson_rd(x, y, cplx<T>(T(1)));
    if (m != T()) r += T(2) * m * comp_ellint_2k(cplx<T>(k));
    return r;
}
template <class T>
[[nodiscard]] inline cplx<T> ellint_3n(T k, const cplx<T>& nu, T phi) {
    const T m = std::round(phi / pi<T>), sg = period_sign(m);
    const T s = sg * detail::lm::sin(phi), c = sg * detail::lm::cos(phi);
    const cplx<T> x(c * c), y(ellint_delta2(k, c));
    cplx<T> r = s * carlson_rf(x, y, cplx<T>(T(1))) +
                (nu / T(3)) * (s * s * s) * carlson_rj(x, y, cplx<T>(T(1)), rj_param_sin(nu, s, c));
    if (m != T()) r += T(2) * m * comp_ellint_3n(k, nu);
    return r;
}

// ---- exponential integrals ----------------------------------------------------------------
// E1(z), principal branch (cut along the negative real axis).  Power series where the
// cancellation exp(|z| + Re z) is small, otherwise the continued fraction
//   E1(z) = e^{-z} / (z + 1 - 1/(z + 3 - 4/(z + 5 - ...)))   (valid |ph z| < pi).
template <class T>
[[nodiscard]] inline cplx<T> expint_e1(const cplx<T>& z) {
    if (!is_finite(z)) {
        if (z.imag() == T() && z.real() == std::numeric_limits<T>::infinity()) return {};
        return {nan<T>(), nan<T>()};
    }
    if (z == cplx<T>()) return {std::numeric_limits<T>::infinity(), T()};
    const T r = std::abs(z);
    if (r <= T(2) || r + z.real() <= T(4)) {
        cplx<T> term(T(1)), sum{};
        for (int k = 1; k < 100000; ++k) {
            term *= -z / T(k);
            const cplx<T> t = term / T(k);
            sum += t;
            if (T(k) > r && std::abs(t) <= eps<T>() * std::abs(sum) / 4) break;
        }
        return -egamma_v<T> - std::log(z) - sum;
    }
    const T tiny = std::numeric_limits<T>::min() * T(1e10);
    cplx<T> b = z + T(1), c = T(1) / tiny, d = cdiv(T(1), b), h = d;
    for (int i = 1; i < 1000000; ++i) {
        const T an = -T(i) * T(i);
        b += T(2);
        d = an * d + b;
        if (std::abs(d) < tiny) d = tiny;
        d = cdiv(T(1), d);
        c = b + cdiv(cplx<T>(an), c);
        if (std::abs(c) < tiny) c = tiny;
        const cplx<T> del = c * d;
        h *= del;
        if (std::abs(del - T(1)) <= eps<T>()) break;
    }
    return h * detail::lm::exp(-z);
}
// Ei(z): on the real axis the real-valued (principal value) Ei of <cmath>; elsewhere
// Ei(z) = -E1(-z) + i pi sgn(Im z), i.e. gamma + log z + sum z^k/(k k!).
template <class T>
[[nodiscard]] inline cplx<T> expint_ei(const cplx<T>& z) {
    if (z.imag() == T()) {
        const T x = z.real();
        if (x == T()) return {-std::numeric_limits<T>::infinity(), T()};
        if (x < T()) return {-expint_e1(cplx<T>(-x)).real(), T()};
        return {-expint_e1(cplx<T>(-x, -T())).real(), T()};
    }
    const cplx<T> e = -expint_e1(-z);
    return {e.real(), e.imag() + (z.imag() > T() ? pi<T> : -pi<T>)};
}

// Si and Ci (Ci has its cut along the negative real axis).
template <class T>
[[nodiscard]] inline cplx<T> sici(const cplx<T>& z, bool want_ci) {
    if (!is_finite(z)) return {nan<T>(), nan<T>()};
    if (z == cplx<T>()) return want_ci ? cplx<T>(-std::numeric_limits<T>::infinity(), T()) : cplx<T>();
    const T r = std::abs(z);
    if (r <= T(4)) {
        const cplx<T> z2 = z * z;
        if (!want_ci) {
            cplx<T> term = z, sum = z;
            for (int k = 1; k < 1000; ++k) {
                term *= -z2 / (T(2 * k) * T(2 * k + 1));
                const cplx<T> t = term / T(2 * k + 1);
                sum += t;
                if (std::abs(t) <= eps<T>() * std::abs(sum) / 4) break;
            }
            return sum;
        }
        cplx<T> term(T(1)), sum{};
        for (int k = 1; k < 1000; ++k) {
            term *= -z2 / (T(2 * k - 1) * T(2 * k));
            const cplx<T> t = term / T(2 * k);
            sum += t;
            if (std::abs(t) <= eps<T>() * std::max(std::abs(sum), T(1)) / 4) break;
        }
        return egamma_v<T> + std::log(z) + sum;
    }
    if (z.real() < T()) {
        // Si(-z) = -Si(z); Ci(z) = Ci(-z) +- i pi for z in the upper/lower half plane.
        const cplx<T> v = sici(-z, want_ci);
        if (!want_ci) return -v;
        return {v.real(), v.imag() + (on_lower_side(z) ? -pi<T> : pi<T>)};
    }
    // |ph z| <= pi/2: DLMF 6.5.5-6.5.6 (signed zeros select the side at ph z = +-pi/2).
    const cplx<T> ep = expint_e1(mul_i(z)), em = expint_e1(mul_mi(z));
    if (want_ci) return -(ep + em) / T(2);
    return pi<T> / T(2) + mul_mi(ep - em) / T(2);
}

// Fresnel integrals S(z) = int_0^z sin(pi t^2/2) dt, C(z) = int_0^z cos(pi t^2/2) dt.
template <class T>
[[nodiscard]] inline cplx<T> fresnel(const cplx<T>& z, bool want_s) {
    if (!is_finite(z)) return {nan<T>(), nan<T>()};
    if (std::abs(z) <= T(1.5)) {
        const cplx<T> z4 = z * z * z * z;
        const T c = pi<T> / T(2);
        // S: sum (-1)^k c^{2k+1} z^{4k+3}/((2k+1)! (4k+3)); C: sum (-1)^k c^{2k} z^{4k+1}/((2k)! (4k+1))
        cplx<T> p = want_s ? c * z * z * z : z, sum = want_s ? p / T(3) : p;
        for (int k = 1; k < 500; ++k) {
            const T a = want_s ? T(2 * k) * T(2 * k + 1) : T(2 * k - 1) * T(2 * k);
            p *= -(c * c) * z4 / a;
            const cplx<T> t = p / T(want_s ? 4 * k + 3 : 4 * k + 1);
            sum += t;
            if (std::abs(t) <= eps<T>() * std::abs(sum) / 4) break;
        }
        return sum;
    }
    // C + iS = (1+i)/2 erf((1-i) sqrt(pi) z/2)   (DLMF 7.5.?) written for C and S separately.
    const T h = detail::lm::sqrt(pi<T>) / T(2);
    const cplx<T> ea = erf(cplx<T>(h, h) * z), eb = erf(cplx<T>(h, -h) * z);
    if (want_s) return cplx<T>(T(0.25), T(0.25)) * (ea - mul_i(eb));
    return cplx<T>(T(0.25), -T(0.25)) * (ea + mul_i(eb));
}

// ---- polygamma ---------------------------------------------------------------------------
template <class T>
[[nodiscard]] inline cplx<T> polygamma_k(unsigned m, cplx<T> z) {
    if (m == 0) return digamma(z);
    if (!is_finite(z)) return {nan<T>(), nan<T>()};
    if (is_nonpositive_integer(z)) return {nan<T>(), T()};
    const T sign = (m % 2 == 1) ? T(1) : T(-1);  // (-1)^{m+1}
    T mfact = T(1);
    for (unsigned i = 2; i <= m; ++i) mfact *= T(i);
    const T R = T(25) + T(m);
    cplx<T> s{};
    while (z.real() < R) {
        s += cdiv(T(1), ipow(z, m + 1));
        z += T(1);
    }
    // psi^(m)(w) ~ (-1)^{m+1} [ (m-1)!/w^m + m!/(2 w^{m+1}) + sum_k B_2k (2k+m-1)!/((2k)! w^{2k+m}) ]
    const cplx<T> inv = cdiv(T(1), z), inv2 = inv * inv;
    cplx<T> wm = ipow(inv, m);  // w^{-m}
    cplx<T> asym = (mfact / T(m)) * wm + (mfact / T(2)) * wm * inv;
    T coef = mfact / T(m);  // (2k+m-1)!/(2k)! for k = 0 is (m-1)!
    cplx<T> p = wm;
    for (int k = 1; k <= stirling_terms; ++k) {
        const T kk = T(k), mm = T(m);
        coef *= (T(2) * kk + mm - T(2)) * (T(2) * kk + mm - T(1)) / ((T(2) * kk - T(1)) * (T(2) * kk));
        p *= inv2;
        asym += bernoulli_2k<T>(k) * coef * p;
    }
    return sign * (asym + mfact * s);
}

// ---- Hurwitz zeta by Hermite's integral (any s != 1, Re a > 0) ------------------------
//   zeta(s, a) = a^{-s}/2 + a^{1-s}/(s-1) + 2 int_0^inf sin(s atan(t/a)) (a^2+t^2)^{-s/2} / (e^{2 pi t} - 1) dt
// (DLMF 25.11.29), with sin(s atan(t/a)) (a^2+t^2)^{-s/2} = ((a - it)^{-s} - (a + it)^{-s}) / (2i).
// No partial sums are formed, so it does not suffer the cancellation that Euler-Maclaurin
// has for Re s < 0.
template <class T>
[[nodiscard]] inline cplx<T> hurwitz_zeta_hermite(const cplx<T>& s, cplx<T> a) {
    // Bring Re a into (lo - 1, lo]: for small |Im s| (lo = 1) the terms below are then of
    // the size of the result.
    cplx<T> head{};
    while (!(a.real() > T())) {
        if (a == cplx<T>()) return {std::numeric_limits<T>::infinity(), T()};
        head += detail::lm::exp(-s * std::log(a));
        a += T(1);
    }
    // For |Im s| large the integrand grows like e^{|Im s| atan(t/a)}; with
    // Re a >= 1.3 |Im s|/(2 pi) that growth is dominated by e^{-2 pi t}, so shift a up to there.
    const T lo = std::max(T(1), std::ceil(T(1.3) * std::abs(s.imag()) / (T(2) * pi<T>)));
    while (a.real() > lo) {
        a -= T(1);
        head -= detail::lm::exp(-s * std::log(a));
    }
    while (a.real() <= lo - T(1)) {
        head += detail::lm::exp(-s * std::log(a));
        a += T(1);
    }
    const T twopi = T(2) * pi<T>;
    const cplx<T> I = de_half_infinite<T>([&](T t) -> cplx<T> {
        // (a - it)^{-s} - (a + it)^{-s} = e^{-sc} (e^{sd} - e^{-sd}), with
        //   c = (log(a - it) + log(a + it))/2,  d = (log(a + it) - log(a - it))/2 = atanh(it/a),
        // d from log1p(w) - log1p(-w) (no cancellation: the two terms have opposite signs).
        // Subtracting the two powers directly cancels like t -> 0, and is exact only if the
        // two products s*log(a -+ it) round symmetrically, which FMA contraction breaks
        // (a 1-ulp asymmetry divided by e^{2 pi t} - 1 ~ 2 pi t overflowed to -inf).
        // Real a (the usual case): c = log|a + it|, d = i atan(t/a).
        cplx<T> c, d;
        if (a.imag() == T()) {
            c = std::log(std::hypot(a.real(), t));
            d = {T(), std::atan2(t, a.real())};
        } else {
            const cplx<T> w = cdiv(cplx<T>(T(), t), a);
            c = (std::log(cplx<T>(a.real(), a.imag() - t)) + std::log(cplx<T>(a.real(), a.imag() + t))) / T(2);
            d = (log1p(w) - log1p(-w)) / T(2);
        }
        // for larger t divide by e^{2 pi t} - 1 inside the exponentials, so that the
        // factors e^{|Im s| atan(t/a)} cannot overflow
        const bool big = t > T(1);
        const T ld = big ? twopi * t + std::log1p(-detail::lm::exp(-twopi * t)) : T();
        const cplx<T> X = -s * c - ld, Y = s * d;  // (c real or not, d imaginary or not)
        // e^{X+Y} - e^{X-Y}: as 2 e^X sinh(Y) where that difference is small (|Re Y| <= 1);
        // otherwise the two terms differ in size and are formed separately (no overflow).
        const cplx<T> diff = std::abs(Y.real()) <= T(1) ? T(2) * detail::lm::exp(X) * detail::lm::sinh(Y)
                                                        : detail::lm::exp(X + Y) - detail::lm::exp(X - Y);
        return big ? mul_mi(diff) : mul_mi(diff) / detail::lm::expm1(twopi * t);
    });
    const cplx<T> la = std::log(a);
    return head + detail::lm::exp(-s * la) / T(2) + cdiv(detail::lm::exp((T(1) - s) * la), s - T(1)) + I;
}

// ---- Hurwitz zeta (complex s, complex a) by Euler-Maclaurin -------------------------------
template <class T>
[[nodiscard]] inline cplx<T> hurwitz_zeta_em(const cplx<T>& s, const cplx<T>& a) {
    if (!is_finite(s) || !is_finite(a)) return {nan<T>(), nan<T>()};
    if (s == cplx<T>(T(1))) return {std::numeric_limits<T>::infinity(), T()};
    const int M = stirling_terms;
    const T need = (std::abs(s) + T(2 * M)) / pi<T> + T(2);
    long N = 0;
    while (std::abs(a + T(N)) < need || (a.real() + T(N)) < T(1)) ++N;
    cplx<T> sum{};
    for (long k = 0; k < N; ++k) {
        const cplx<T> ak = a + T(k);
        if (ak == cplx<T>()) return {std::numeric_limits<T>::infinity(), T()};
        sum += detail::lm::exp(-s * std::log(ak));
    }
    const cplx<T> w = a + T(N), lw = std::log(w);
    const cplx<T> ws = detail::lm::exp(-s * lw);  // w^{-s}
    cplx<T> tail = cdiv(w * ws, s - T(1)) + ws / T(2);
    // sum_j B_2j/(2j)! (s)_{2j-1} w^{-s-2j+1}
    const cplx<T> inv = cdiv(T(1), w), inv2 = inv * inv;
    cplx<T> poch = s, p = ws * inv;  // (s)_1 w^{-s-1}
    T fact = T(2);                    // (2j)!
    for (int j = 1; j <= M; ++j) {
        tail += bernoulli_2k<T>(j) / fact * poch * p;
        poch *= (s + T(2 * j - 1)) * (s + T(2 * j));
        p *= inv2;
        fact *= T(2 * j + 1) * T(2 * j + 2);
    }
    return sum + tail;
}

// ---- Hurwitz zeta by Hurwitz's formula (real a, Re s <= -2) -------------------------------
//   zeta(s, a) = Gamma(t)/(2 pi)^t [e^{i pi t/2} sum e^{-2 pi i n a} n^{-t} + e^{-i pi t/2} sum e^{2 pi i n a} n^{-t}],
// t = 1 - s, 0 < a <= 1 (DLMF 25.11.9).  The two terms do not cancel, so this stays accurate
// for large |Im s| where Hermite's integral cancels.  The sums converge like n^{-Re t}.
template <class T>
[[nodiscard]] inline cplx<T> hurwitz_zeta_reflect(const cplx<T>& s, T a) {
    cplx<T> head{};
    while (!(a > T())) {
        head += detail::lm::exp(-s * std::log(cplx<T>(a)));
        a += T(1);
    }
    while (a > T(1)) {
        a -= T(1);
        head -= detail::lm::exp(-s * std::log(cplx<T>(a)));
    }
    const cplx<T> t = T(1) - s;
    cplx<T> sp{}, sm{};  // sum e^{+2 pi i n a} n^{-t}, sum e^{-2 pi i n a} n^{-t}
    const T stop = eps<T>() / T(16);
    for (long n = 1; n < 100000000; ++n) {
        const T ln = std::log(T(n));
        const cplx<T> nt = detail::lm::exp(-t * ln);
        const T x = T(2) * T(n) * a;  // phase / pi, reduced exactly-ish via sinpi/cospi
        const T fr = x - T(2) * std::floor(x / T(2));
        const cplx<T> e = expi_pi(fr);
        sp += e * nt;
        sm += std::conj(e) * nt;
        if (detail::lm::exp(-t.real() * ln) <= stop * std::max(std::abs(sp), std::abs(sm)) && n > 8) break;
    }
    const cplx<T> lg = lgamma(t) - t * std::log(T(2) * pi<T>);
    const cplx<T> ip(T(), pi<T> / T(2));
    return head + detail::lm::exp(lg + ip * t) * sm + detail::lm::exp(lg - ip * t) * sp;
}

template <class T>
[[nodiscard]] inline cplx<T> hurwitz_zeta_k(const cplx<T>& s, const cplx<T>& a) {
    if (!is_finite(s) || !is_finite(a)) return {nan<T>(), nan<T>()};
    if (s == cplx<T>(T(1))) return {std::numeric_limits<T>::infinity(), T()};
    // Euler-Maclaurin's partial sums cancel for Re s < 1/2, and its fixed number of Bernoulli
    // terms limits it for large |Im s|.  Hurwitz's formula covers real a with Re s <= -2 and
    // |Im s| > 8.  Hermite's integral (with a shifted so the integrand does not grow with
    // |Im s|) covers the rest; it cancels only for Re s > 1/2 with small |Im s|, where
    // Euler-Maclaurin is used instead.
    const T ai = std::abs(s.imag());
    if (ai > T(8) && s.real() <= T(-2) && a.imag() == T()) return hurwitz_zeta_reflect(s, a.real());
    if (s.real() >= T(0.5) && ai <= (a.imag() == T() ? T(15) : T(4))) return hurwitz_zeta_em(s, a);
    return hurwitz_zeta_hermite(s, a);
}

// ---- polylogarithm Li_s(z), real s ----------------------------------------------------------
// Bernoulli polynomial B_n(x), n <= 40.
template <class T>
[[nodiscard]] inline cplx<T> bernoulli_poly(unsigned n, const cplx<T>& x) {
    // B_n(x) = sum_k C(n,k) B_k x^{n-k}, B_0 = 1, B_1 = -1/2, B_odd>1 = 0
    cplx<T> sum{}, xp(T(1));  // xp = x^{n-k} built from k = n downwards
    T binom = T(1);           // C(n, k) for k = n
    for (unsigned k = n + 1; k-- > 0;) {
        T bk = T();
        if (k == 0) bk = T(1);
        else if (k == 1) bk = T(-0.5);
        else if (k % 2 == 0) bk = bernoulli_2k<T>(static_cast<int>(k / 2));
        sum += bk * binom * xp;
        xp *= x;
        if (k > 0) binom = binom * T(k) / T(n - k + 1);  // C(n, k-1)
    }
    return sum;
}
template <class T>
[[nodiscard]] inline cplx<T> polylog_k(T s, const cplx<T>& z);

template <class T>
[[nodiscard]] inline cplx<T> polylog_series(T s, const cplx<T>& z) {
    cplx<T> p = z, sum{};
    for (int k = 1; k < 100000; ++k) {
        const cplx<T> t = p * detail::lm::pow(T(k), -s);
        sum += t;
        if (std::abs(t) <= eps<T>() * std::abs(sum) / 8) break;
        p *= z;
    }
    return sum;
}
// Expansion in mu = log z (|mu| < 2 pi):
//   non-integer s: Li_s = Gamma(1-s) (-mu)^{s-1} + sum_k zeta(s-k) mu^k/k!
//   integer s = n >= 1: Li_n = mu^{n-1}/(n-1)! (H_{n-1} - log(-mu)) + sum_{k != n-1} zeta(n-k) mu^k/k!
template <class T>
[[nodiscard]] inline cplx<T> polylog_logseries(T s, const cplx<T>& z) {
    const cplx<T> mu = std::log(z);
    const bool pos_int = is_integer(s) && s >= T(1);
    // s = n + delta next to a positive integer n: the k = n-1 term zeta(1+delta) mu^{n-1}/(n-1)!
    // and Gamma(1-s)(-mu)^{s-1} both have a 1/delta pole.  They are combined analytically:
    //   mu^{n-1}/(n-1)! [A(delta) - expm1(E)/delta],  A = zeta(1+delta) - 1/delta (Stieltjes),
    //   E = delta log(-mu) + log(pi delta/sin(pi delta)) + log Gamma(n) - log Gamma(n+delta),
    // with E/delta from its Taylor series, so nothing cancels.
    const T sr = std::nearbyint(s), delta = s - sr;
    const bool near_int = !pos_int && sr >= T(1) && std::abs(delta) <= T(0.1);
    const long n = pos_int ? static_cast<long>(s) : near_int ? static_cast<long>(sr) : 0;
    cplx<T> sum{}, p(T(1));
    T fact = T(1);
    for (long k = 0; k < 2000; ++k) {
        if (k > 0) {
            p *= mu;
            fact *= T(k);
        }
        if ((pos_int || near_int) && k == n - 1) continue;
        const T sk = s - T(k);
        if (sk < T() && is_integer(sk) && static_cast<long>(sk) % 2 == 0) continue;  // trivial zeros
        // Next to a trivial zero -2m, zeta(sk) ~ eta zeta'(-2m) with eta = sk + 2m; the rounding
        // of s - k (ulp(k), far above eta for s just off an integer) must not enter eta.  So
        // zeta(sk) = (-1)^m sin(pi eta/2) 2^sk pi^(sk-1) Gamma(1-sk) zeta(1-sk), eta exact.
        const T sh = sk / T(2) - std::nearbyint(sk / T(2));
        const bool near_trivial = sk < T(-1) && std::abs(sh) < T(0.05);
        T zk;
        if (near_trivial) {
            const T m = -std::nearbyint(sk / T(2));
            const T eta = ((sr - T(k)) + T(2) * m) + (s - sr);  // exact
            const T mag = detail::lm::exp(sk * std::log(T(2)) + (sk - T(1)) * std::log(pi<T>) + detail::lm::lgamma(T(1) - sk)) *
                          riemann_zeta(cplx<T>(T(1) - sk)).real();
            zk = (std::fmod(m, T(2)) == T() ? T(1) : T(-1)) * sinpi(eta / T(2)) * mag;
        } else {
            zk = riemann_zeta(cplx<T>(sk)).real();
        }
        const cplx<T> t = zk * p / fact;
        sum += t;
        // A term next to a trivial zero (s - k near -2, -4, ...) is small without the series
        // having converged: do not stop there (this cut Li_s short for s near an integer).
        if (k > 4 && !near_trivial && std::abs(t) <= eps<T>() * std::abs(sum) / 8) break;
    }
    if (pos_int) {
        T H = T();
        for (long i = 1; i < n; ++i) H += T(1) / T(i);
        T f = T(1);
        for (long i = 2; i < n; ++i) f *= T(i);
        return sum + ipow(mu, static_cast<unsigned>(n - 1)) / f * (H - std::log(-mu));
    }
    if (near_int) {
        static constexpr long double stieltjes_a[] = {  // (-1)^k gamma_k / k!
        0.577215664901532860606512090082402431L,
        0.0728158454836767248605863758749013191L,
        -0.00484518159643615924226519301760626468L,
        -0.000342305736717224311026674423792230714L,
        0.0000968904193944708357278404240635861667L,
        -0.00000661103181084218918127779064537037359L,
        -0.000000331624090875277235933919224872228164L,
        0.000000104620945844791874221051092356916326L,
        -0.00000000873321810027379736116200806054585173L,
        0.00000000009478277782762358955554074018470952L,
        0.0000000000565842192760870796637241758259108529L,
        -6.7686898635136966558667549018948331e-12L,
        3.49211593667203185445521566819783998e-13L,
        4.41042474175775338023724020535713989e-15L,
        -2.39978622177099917550505620996371351e-15L,
        2.16773122007268285496389323582067013e-16L,
        -9.54446607636696517342499076972236294e-18L,
        -7.3876766605386364978155782870704015e-20L,
        4.80085078248806522761766273828542567e-20L,
        -4.13995673771330564126948200778884778e-21L,
        1.91682015939912339496481976180313116e-22L,
        -2.0441543122262166077275863759982535e-24L,
        -4.81849850110735344392921837431093251e-25L,
        4.81185705151256647946110663519045342e-26L,
        -2.56026331031881493660912995816108722e-27L,
        6.9278408953046671238801264492378066e-29L,
        1.62860755048558674407104276822745755e-30L,
        -3.19393756115325557604210739567211501e-31L,
        2.09915158936342552768548584185771825e-32L,
        -8.33674529544144047562508226550865547e-34L
};
        T A{}, dk = T(1);
        for (const long double a : stieltjes_a) {
            A += static_cast<T>(a) * dk;
            dk *= delta;
        }
        // E/delta = log(-mu) - psi(n) + sum_{k>=1} zeta(2k) delta^{2k-1}/k
        //           - sum_{j>=2} (-1)^j (zeta(j) - H^{(j)}_{n-1}) delta^{j-1}/j
        T psi = -static_cast<T>(stieltjes_a[0]);
        for (long i = 1; i < n; ++i) psi += T(1) / T(i);
        T ser{}, dj = T(1);  // dj = delta^{j-1}
        for (int j = 2; j < 60; ++j) {
            dj *= delta;
            const T zj = riemann_zeta(cplx<T>(T(j))).real();
            T Hj{};
            for (long i = 1; i < n; ++i) Hj += detail::lm::pow(T(i), -T(j));
            T t = -((j % 2 == 0) ? T(1) : T(-1)) * (zj - Hj) / T(j);
            if (j % 2 == 0) t += zj / T(j / 2);  // zeta(2k)/k, k = j/2
            t *= dj;
            ser += t;
            if (std::abs(dj) * T(2) <= eps<T>() / T(16)) break;  // coefficients are below 2
        }
        const cplx<T> Ed = std::log(-mu) - psi + ser;
        const cplx<T> E = delta * Ed;
        // expm1(E)/E, accurate for small E
        cplx<T> q;
        if (std::abs(E) < T(1e-3)) {
            cplx<T> term(T(1)), acc(T(1));
            for (int i = 2; i < 30; ++i) {
                term *= E / T(i);
                acc += term;
                if (std::abs(term) <= eps<T>() / T(8)) break;
            }
            q = acc;
        } else {
            const T c = detail::lm::cos(E.imag());
            const cplx<T> em1(detail::lm::expm1(E.real()) * c - T(2) * detail::lm::pow(detail::lm::sin(E.imag() / T(2)), 2),
                              detail::lm::exp(E.real()) * detail::lm::sin(E.imag()));
            q = cdiv(em1, E);
        }
        T f = T(1);
        for (long i = 2; i < n; ++i) f *= T(i);
        return sum + ipow(mu, static_cast<unsigned>(n - 1)) / f * (A - Ed * q);
    }
    return sum + detail::lm::exp(lgamma(cplx<T>(T(1) - s))) * detail::lm::exp((s - T(1)) * std::log(-mu));
}
template <class T>
[[nodiscard]] inline cplx<T> polylog_k(T s, const cplx<T>& z) {
    if (!is_finite(z) || !std::isfinite(s)) return {nan<T>(), nan<T>()};
    if (z == cplx<T>()) return {};
    if (z == cplx<T>(T(1))) {
        if (s > T(1)) return riemann_zeta(cplx<T>(s));
        return {std::numeric_limits<T>::infinity(), T()};
    }
    if (s == T()) return cdiv(z, T(1) - z);
    const T r = std::abs(z);
    if (r <= T(0.5)) return polylog_series(s, z);
    if (r >= T(2)) {
        const cplx<T> zi = cdiv(T(1), z);
        if (is_integer(s) && s < T()) {  // Li_{-n}(z) = (-1)^{n+1} Li_{-n}(1/z)
            const long nn = static_cast<long>(-s);
            const cplx<T> v = polylog_series(s, zi);
            return (nn % 2 == 0) ? -v : v;
        }
        // a = 1/2 + log(-z)/(2 pi i)
        const cplx<T> a = T(0.5) + mul_mi(std::log(-z)) / (T(2) * pi<T>);
        if (is_integer(s) && s <= T(40)) {
            // Li_n(z) = -(-1)^n Li_n(1/z) - (2 pi i)^n/n! B_n(a)
            const unsigned nn = static_cast<unsigned>(s);
            T f = T(1);
            for (unsigned i = 2; i <= nn; ++i) f *= T(i);
            const cplx<T> c = ipow(cplx<T>(T(), T(2) * pi<T>), nn) / f;
            const cplx<T> v = polylog_series(s, zi);
            return ((nn % 2 == 0) ? -v : v) - c * bernoulli_poly(nn, a);
        }
        if (!is_integer(s) && (s <= T(1) || std::abs(std::log(z)) > T(5.5))) {
            // Li_s(z) = -e^{i pi s} Li_s(1/z) + (2 pi)^s e^{i pi s/2}/Gamma(s) zeta(1-s, a)
            const cplx<T> v = polylog_series(s, zi);
            const cplx<T> g = detail::lm::exp(s * std::log(T(2) * pi<T>) - lgamma(cplx<T>(s))) * expi_pi(s / T(2));
            return -expi_pi(s) * v + g * hurwitz_zeta_k(cplx<T>(T(1) - s), a);
        }
    }
    return polylog_logseries(s, z);
}

// 1/e as hi + lo (each rounded to T), for evaluating z + 1/e near the Lambert W branch point.
template <class T>
[[nodiscard]] inline T inv_e_hi() noexcept {
    if constexpr (std::numeric_limits<T>::digits <= 53) return static_cast<T>(0x1.78b56362cef38p-2);
    else return static_cast<T>(0xbc5ab1b16779be35p-65L);
}
template <class T>
[[nodiscard]] inline T inv_e_lo() noexcept {
    if constexpr (std::numeric_limits<T>::digits <= 53) return static_cast<T>(-0x1.ca8a4270fadf5p-57);
    else return static_cast<T>(0xeb7b1e0a4153e437p-130L);
}
// g(d) = (d - 1) e^d + 1 = sum_{n>=2} (n-1) d^n / n!, accurate for small d.
template <class T>
[[nodiscard]] inline cplx<T> lambert_g(const cplx<T>& d) {
    if (std::abs(d) > T(0.5)) return (d - T(1)) * detail::lm::exp(d) + T(1);
    cplx<T> term = d * d / T(2), sum = term;  // n = 2: (n-1)/n! = 1/2
    for (int n = 3; n < 60; ++n) {
        term *= d / T(n);  // d^n / n!
        const cplx<T> t = term * T(n - 1);
        sum += t;
        if (std::abs(t) <= eps<T>() * std::abs(sum) / 8) break;
    }
    return sum;
}

// ---- Lambert W ------------------------------------------------------------------------------
template <class T>
[[nodiscard]] inline cplx<T> lambert_w_k(int k, const cplx<T>& z) {
    if (!is_finite(z)) return {nan<T>(), nan<T>()};
    if (z == cplx<T>()) {
        if (k == 0) return {};
        return {-std::numeric_limits<T>::infinity(), T()};
    }
    const T e1 = detail::lm::exp(T(-1));
    const bool real_axis = z.imag() == T();
    cplx<T> w;
    // Distance from the branch point, z + 1/e, with 1/e split into hi + lo so the sum is
    // exact up to one final rounding (z + hi is exact near -1/e by Sterbenz's lemma).
    const cplx<T> zb = (z + inv_e_hi<T>()) + inv_e_lo<T>();
    const bool near_bp = std::abs(zb) < T(0.25) || (real_axis && z.real() < -e1 && z.real() > T(-1));
    // p = sqrt(2(e z + 1)); W = -1 + p - p^2/3 + 11 p^3/72 - ... on the branch selected by p's sign.
    auto bp_series = [&](cplx<T> p) {
        return T(-1) + p - p * p / T(3) + T(11) / T(72) * p * p * p - T(43) / T(540) * p * p * p * p;
    };
    const cplx<T> p = detail::lm::sqrt(T(2) * detail::lm::exp(T(1)) * zb);
    if (near_bp && k == 0) {
        w = bp_series(p);
    } else if (near_bp && k == -1 && !on_lower_side(z)) {
        w = bp_series(-p);
    } else if (near_bp && k == 1 && on_lower_side(z)) {
        w = bp_series(-p);
    } else if (k == 0 && std::abs(z) <= T(0.3)) {
        w = z * (T(1) - z + T(1.5) * z * z);  // series start
    } else if (k == 0 && std::abs(z) <= T(3)) {
        // log(1 + z), or the branch-point series with the principal root near z = -1.  (The
        // asymptotic start log z - log log z used to be taken for 1 < |z| <= 3, where log z is
        // small and Halley's iteration could converge to W_{+-1}, e.g. at z = 0.97 + 0.27i.)
        w = std::abs(T(1) + z) < T(0.5) ? bp_series(p) : std::log(T(1) + z);
    } else {
        cplx<T> L1 = std::log(z);
        L1 += cplx<T>(T(), T(2) * pi<T> * T(k));
        if (real_axis && k == -1 && z.real() < T() && z.real() > -e1) L1 = std::log(-z) - std::log(-std::log(-z));
        const cplx<T> L2 = std::log(L1);
        w = (std::abs(L1) > T(0)) ? L1 - L2 + cdiv(L2, L1) : L1;
        if (real_axis && k == -1 && z.real() < T() && z.real() > -e1) w = L1;  // real start
    }
    auto halley = [&](cplx<T> x) {
        for (int it = 0; it < 200; ++it) {
            cplx<T> dw;
            if (std::abs(x) > T(30) && std::abs(z) > T(1e6)) {
                // Newton on w + log w = log z + 2 pi i k (avoids overflow of w e^w).
                const cplx<T> g = x + std::log(x) - std::log(z) - cplx<T>(T(), T(2) * pi<T> * T(k));
                dw = cdiv(g, T(1) + cdiv(T(1), x));
            } else {
                const cplx<T> ew = detail::lm::exp(x), f = x * ew - z, w1 = x + T(1);
                const cplx<T> den = ew * w1 - cdiv((x + T(2)) * f, T(2) * w1);
                dw = cdiv(f, den);
            }
            if (!is_finite(dw)) break;  // at the branch point w = -1 exactly
            x -= dw;
            if (!(std::abs(dw) > T(4) * eps<T>() * std::abs(x))) break;
        }
        return x;
    };
    // W_0 takes values left-bounded by the curve -eta cot(eta) + i eta, |eta| < pi.
    auto in_w0 = [](const cplx<T>& v) {
        const T eta = v.imag();
        if (!(std::abs(eta) < pi<T>)) return false;
        if (eta == T()) return v.real() >= T(-1) - T(1e-6);
        return v.real() > -eta / detail::lm::tan(eta) - T(1e-6) * (T(1) + std::abs(v.real()));
    };
    w = halley(w);
    if (k == 0 && !near_bp && !in_w0(w)) {
        // The start led Halley's iteration to another branch: retry from other starts.
        const cplx<T> lz = std::log(T(1) + z);
        const cplx<T> starts[] = {bp_series(p), lz * (T(1) - cdiv(std::log(T(1) + lz), T(2) + lz)), lz, cplx<T>(T(0.5))};
        for (const auto& s0 : starts) {
            const cplx<T> v = halley(s0);
            if (in_w0(v)) {
                w = v;
                break;
            }
        }
    }
    // Near the branch point w e^w - z cancels badly.  Refine d = w + 1 by Newton on the
    // well-conditioned residual  e^{-1} g(d) - (z + 1/e),  g(d) = (d-1) e^d + 1.
    cplx<T> d = w + T(1);
    if (is_finite(d) && std::abs(d) < T(0.5)) {
        for (int it = 0; it < 8; ++it) {
            const cplx<T> r = lambert_g(d) * inv_e_hi<T>() + lambert_g(d) * inv_e_lo<T>() - zb;
            const cplx<T> der = detail::lm::exp(d) * d * e1;  // d/dd [e^{-1} g(d)] = e^{d-1} d
            if (der == cplx<T>()) break;
            const cplx<T> step = cdiv(r, der);
            if (!is_finite(step)) break;
            d -= step;
            if (!(std::abs(step) > eps<T>() * std::abs(d))) break;
        }
        w = d - T(1);
    }
    return w;
}

// ---- incomplete gamma (real a, complex z) ---------------------------------------------------
// Returns {P, Q} regularized, for a > 0.
// exp(a log z + c z) for real a and c = 0 or -1, with the exponent in double-T so that
// large |a log z| or |z| do not amplify its rounding.
template <class T>
[[nodiscard]] inline cplx<T> exp_alogz(T a, const cplx<T>& z, T c) noexcept {
    const cdd<T> lz = clog_dd(dd<T>{z.real(), T()}, dd<T>{z.imag(), T()});
    cdd<T> e{lz.re * a, lz.im * a};
    if (c != T()) {
        e.re = e.re + z.real() * c;
        e.im = e.im + z.imag() * c;
    }
    return exp_split(e);
}
template <class T>
struct pq { cplx<T> p, q; };
template <class T>
[[nodiscard]] inline pq<T> gamma_pq(T a, const cplx<T>& z) {
    if (z == cplx<T>()) return {cplx<T>(), cplx<T>(T(1))};
    const T r = std::abs(z);
    const cplx<T> lz = std::log(z);
    const T lga = detail::lm::lgamma(a);
    if (r < a + T(1)) {
        // P = z^a e^{-z}/Gamma(a) sum z^k / (a (a+1) ... (a+k))
        cplx<T> term(T(1) / a), sum = term;
        for (int k = 1; k < 1000000; ++k) {
            term *= z / (a + T(k));
            sum += term;
            if (std::abs(term) <= eps<T>() * std::abs(sum) / 4) break;
        }
        const cplx<T> P = detail::lm::exp(a * lz - z - lga) * sum;
        return {P, T(1) - P};
    }
    if (std::abs(std::arg(z)) < T(0.9) * pi<T>) {
        // Q by the Legendre continued fraction.
        const T tiny = std::numeric_limits<T>::min() * T(1e10);
        cplx<T> b = z + T(1) - a, c = T(1) / tiny, d = cdiv(T(1), b), h = d;
        for (int i = 1; i < 1000000; ++i) {
            const T an = -T(i) * (T(i) - a);
            b += T(2);
            d = an * d + b;
            if (std::abs(d) < tiny) d = tiny;
            d = cdiv(T(1), d);
            c = b + cdiv(cplx<T>(an), c);
            if (std::abs(c) < tiny) c = tiny;
            const cplx<T> del = c * d;
            h *= del;
            if (std::abs(del - T(1)) <= eps<T>()) break;
        }
        const cplx<T> Q = detail::lm::exp(a * lz - z - lga) * h;
        return {T(1) - Q, Q};
    }
    // Near the negative real axis: P = z^a/Gamma(a) sum (-z)^k / (k! (a+k)).
    cplx<T> term(T(1)), sum(T(1) / a);
    for (int k = 1; k < 1000000; ++k) {
        term *= -z / T(k);
        const cplx<T> t = term / (a + T(k));
        sum += t;
        if (T(k) > r && std::abs(t) <= eps<T>() * std::abs(sum) / 4) break;
    }
    const cplx<T> P = detail::lm::exp(a * lz - lga) * sum;
    return {P, T(1) - P};
}
// Gamma(a, z) = z^a e^{-z} / (z + 1 - a - 1(1-a)/(z + 3 - a - ...)), |ph z| < pi, any real a.
// The modified Lentz algorithm finds the number of terms N; the fraction is then
// re-evaluated backwards from term N + N/4 + 8, which does not accumulate the forward
// rounding of the Lentz product (several tens of ulp near the negative real axis).
template <class T>
[[nodiscard]] inline cplx<T> gamma_cf(T a, const cplx<T>& z) {
    const T tiny = std::numeric_limits<T>::min() * T(1e10);
    cplx<T> b = z + T(1) - a, c = T(1) / tiny, d = cdiv(T(1), b), h = d;
    int n = 1;
    for (; n < 1000000; ++n) {
        const T an = -T(n) * (T(n) - a);
        b += T(2);
        d = an * d + b;
        if (std::abs(d) < tiny) d = tiny;
        d = cdiv(T(1), d);
        c = b + cdiv(cplx<T>(an), c);
        if (std::abs(c) < tiny) c = tiny;
        const cplx<T> del = c * d;
        h *= del;
        if (std::abs(del - T(1)) <= eps<T>()) break;
    }
    const int m = n + n / 4 + 8;
    const cplx<T> b0 = z + T(1) - a;
    cplx<T> t = b0 + T(2 * m);
    for (int i = m; i >= 1; --i) {
        const T an = -T(i) * (T(i) - a);
        t = b0 + T(2 * (i - 1)) + cdiv(cplx<T>(an), t);
    }
    return cdiv(exp_alogz(a, z, T(-1)), t);
}
// Upper incomplete gamma Gamma(a, z) for any real a (z != 0).
template <class T>
[[nodiscard]] inline cplx<T> gamma_upper(T a, const cplx<T>& z) {
    if (!is_finite(z)) return {nan<T>(), nan<T>()};
    if (z == cplx<T>()) {
        if (a > T()) return {detail::lm::tgamma(a), T()};
        return {std::numeric_limits<T>::infinity(), T()};
    }
    const T az = std::abs(z);
    // The (backward-evaluated) continued fraction is accurate up to about 3 degrees from
    // the negative real axis, i.e. 1 + cos(arg z) > 0.004.  Closer than that, and for
    // |z| <= 2, the alternating power series is used; it cancels only by about
    // e^{|z| + Re z} = e^{0.004 |z|} there.
    // For 0.8 < |z| <= 2 it is also accurate as long as arg z <= 0.6 pi.
    const bool cf_sector = (az > T(2) && az + z.real() > T(0.004) * az) ||
                           (az > T(0.8) && z.real() >= T(-0.3) * az);
    if (a > T(0.5)) {
        if (cf_sector && az >= a + T(1)) return gamma_cf(a, z);
        // Two power series for gamma(a, z), Gamma(a, z) = Gamma(a) - gamma(a, z):
        //   Kummer:      z^a e^{-z} sum z^k / (a (a+1) ... (a+k))
        //   alternating: z^a        sum (-z)^k / (k! (a+k))
        // Each is accepted if its estimated relative error, (largest term / |sum|)
        // |gamma| + Gamma(a), over |Gamma(a, z)|, is at most 16 (in units of eps).  The
        // Kummer form suits Re z >= 0 or |z| < a, the alternating one Re z < 0.
        const T ga = detail::lm::tgamma(a);
        cplx<T> best;
        T best_err = std::numeric_limits<T>::infinity();
        auto consider = [&](const cplx<T>& g, T ratio) {
            const cplx<T> res = ga - g;
            const T err = (ratio * std::abs(g) + ga) / std::abs(res);
            if (!(err >= best_err)) best = res, best_err = err;
        };
        {
            cplx<T> term(T(1) / a), sum = term;
            T mx = std::abs(term);
            for (int k = 1; k < 1000000; ++k) {
                term *= z / (a + T(k));
                sum += term;
                mx = std::max(mx, std::abs(term));
                if (std::abs(term) <= eps<T>() * std::abs(sum) / 4) break;
            }
            consider(exp_alogz(a, z, T(-1)) * sum, mx / std::abs(sum));
        }
        if (best_err <= T(16)) return best;
        {
            cplx<T> term(T(1)), sum(T(1) / a);
            T mx = std::abs(sum);
            for (int k = 1; k < 1000000; ++k) {
                term *= -z / T(k);
                const cplx<T> t = term / (a + T(k));
                sum += t;
                mx = std::max(mx, std::abs(t));
                if (T(k) > az && std::abs(t) <= eps<T>() * std::abs(sum) / 4) break;
            }
            consider(exp_alogz(a, z, T()) * sum, mx / std::abs(sum));
        }
        // The fraction becomes ill-conditioned for |z| well below a; keep the series there.
        if (best_err <= T(16) || !cf_sector || az < T(0.5) * a) return best;
        return gamma_cf(a, z);
    }
    {
        if (cf_sector) return gamma_cf(a, z);
        // Series Gamma(a, z) = Gamma(a) - z^a sum_k (-z)^k / (k! (a+k)), written around the
        // nearest pole a = -n + delta (|delta| <= 1/2).  Gamma(a) and the k = n term
        // both have a pole at delta = 0; combined analytically,
        //   Gamma(a) - z^a (-z)^n / (n! (a+n))
        //     = (-1)^n / n! * [expm1(log f) - expm1(delta log z)] / delta,
        //   f = Gamma(1+delta) / prod_{m=1..n} (1 - delta/m),
        // so no cancellation occurs (the old downward recurrence lost |z|^n-fold).  At
        // delta = 0 (a a non-positive integer) the bracket is H_n - gamma - log z.
        const T nn = std::nearbyint(-a), delta = a + nn;
        const long n = static_cast<long>(nn);
        const cplx<T> lz = std::log(z);
        static constexpr long double zeta_k[] = {
            1.6449340668482264365L, 1.2020569031595942854L, 1.0823232337111381915L,
            1.0369277551433699263L, 1.0173430619844491397L, 1.0083492773819228268L,
            1.0040773561979443394L, 1.0020083928260822144L, 1.0009945751278180853L};
        const T euler_gamma = T(0.57721566490153286060651209L);
        // log Gamma(1 + delta) = -gamma delta + sum_{k>=2} (-1)^k zeta(k) delta^k / k
        // (zeta(k) for k > 10 from its first five terms; |delta| <= 0.1 makes the
        // truncation negligible.)
        T lf = -euler_gamma * delta, dk = -delta;
        if (std::abs(delta) > T(0.1)) {
            lf = detail::lm::lgamma(T(1) + delta);  // the division by delta costs little here
        } else {
            for (int k = 2; k <= 24; ++k) {
                dk *= -delta;
                T zk;
                if (k <= 10) {
                    zk = T(zeta_k[k - 2]);
                } else {
                    zk = T(1);
                    for (int j = 2; j <= 5; ++j) zk += detail::lm::pow(T(j), T(-k));
                }
                lf += zk * dk / T(k);
            }
        }
        T hn = T();
        for (long m = 1; m <= n; ++m) {
            lf -= std::log1p(-delta / T(m));
            hn += T(1) / T(m);
        }
        cplx<T> head;
        if (delta == T()) {
            head = (hn - euler_gamma) - lz;
        } else {
            head = (detail::lm::expm1(lf) - expm1(delta * lz)) / delta;
        }
        T nf = T(1);
        for (long m = 2; m <= n; ++m) nf *= T(m);
        head *= ((n % 2) ? T(-1) : T(1)) / nf;
        // Remaining series terms k != n.
        cplx<T> term(T(1)), sum(n == 0 ? T() : T(1) / a);
        const T r = std::abs(z);
        for (long k = 1; k < 1000000; ++k) {
            term *= -z / T(k);
            if (k == n) continue;
            const cplx<T> t = term / (a + T(k));
            sum += t;
            if (T(k) > r && k > n && std::abs(t) <= eps<T>() * std::abs(sum) / 4) break;
        }
        return head - exp_alogz(a, z, T()) * sum;
    }
}
template <class T>
[[nodiscard]] inline cplx<T> gamma_lower(T a, const cplx<T>& z) {
    if (a <= T() && is_integer(a)) return {nan<T>(), T()};
    if (a > T()) return gamma_pq(a, z).p * detail::lm::tgamma(a);
    return detail::lm::tgamma(a) - gamma_upper(a, z);
}

// ---- regularized incomplete beta I_x(a, b), real a, b > 0, complex x -----------------------
template <class T>
[[nodiscard]] inline cplx<T> betacf(T a, T b, const cplx<T>& x) {
    const T tiny = std::numeric_limits<T>::min() * T(1e10);
    const T qab = a + b, qap = a + T(1), qam = a - T(1);
    cplx<T> c(T(1)), d = T(1) - qab * x / qap;
    if (std::abs(d) < tiny) d = tiny;
    d = cdiv(T(1), d);
    cplx<T> h = d;
    for (int m = 1; m < 1000000; ++m) {
        const T m2 = T(2 * m);
        cplx<T> aa = T(m) * (b - T(m)) * x / ((qam + m2) * (a + m2));
        d = T(1) + aa * d;
        if (std::abs(d) < tiny) d = tiny;
        c = T(1) + cdiv(aa, c);
        if (std::abs(c) < tiny) c = tiny;
        d = cdiv(T(1), d);
        h *= d * c;
        aa = -(a + T(m)) * (qab + T(m)) * x / ((a + m2) * (qap + m2));
        d = T(1) + aa * d;
        if (std::abs(d) < tiny) d = tiny;
        c = T(1) + cdiv(aa, c);
        if (std::abs(c) < tiny) c = tiny;
        d = cdiv(T(1), d);
        const cplx<T> del = d * c;
        h *= del;
        if (std::abs(del - T(1)) <= eps<T>()) break;
    }
    return h;
}
template <class T>
[[nodiscard]] inline cplx<T> ibeta_k(T a, T b, const cplx<T>& x) {
    if (!is_finite(x)) return {nan<T>(), nan<T>()};
    if (x == cplx<T>()) return {};
    if (x == cplx<T>(T(1))) return {T(1), T()};
    const T lbeta = detail::lm::lgamma(a) + detail::lm::lgamma(b) - detail::lm::lgamma(a + b);
    const cplx<T> front = detail::lm::exp(a * std::log(x) + b * log1p(-x) - lbeta);
    if (x.real() < (a + T(1)) / (a + b + T(2))) return front * betacf(a, b, x) / a;
    return T(1) - front * betacf(b, a, T(1) - x) / b;
}

// ---- the axis-reduction lift ----------------------------------------------------------
// F maps std::complex<W> -> std::complex<W>.  If real_on_real is true, the function is
// known to be real on the real axis, and a purely real argument yields a real result.
// f_lower(x) gives f(x - i0), the value below a cut on the real axis; for the functions with
// conjugate symmetry f(conj z) = conj f(z) it is conj f(x + i0).
template <class T, class F, class G>
[[nodiscard]] inline quaternion<T> lift(const quaternion<T>& q, F&& f, bool real_on_real, G&& f_lower) {
    using W = work_t<T>;
    const W v = static_cast<W>(quat::detail::vabs(q));
    if (v == W() && !real_on_real) {
        // On the real axis where f is not real (a cut): the sign of the zero i component
        // chooses the side and the j/k zeros flip when sign(Im f) != sign(i), as in
        // quat::detail::lift_sv for the elementary functions.
        const T b = q.x();
        const cplx<W> r = std::signbit(b) ? f_lower(static_cast<W>(q.w())) : f(cplx<W>(static_cast<W>(q.w()), W()));
        const T im = static_cast<T>(r.imag());
        const bool flip = !std::isnan(im) && std::signbit(im) != std::signbit(b);
        return {static_cast<T>(r.real()), im, flip ? -q.y() : q.y(), flip ? -q.z() : q.z()};
    }
    const cplx<W> r = f(cplx<W>(static_cast<W>(q.w()), v));
    if (v == W()) return {static_cast<T>(r.real()), q.x(), q.y(), q.z()};
    const W s = r.imag() / v;
    // zero components give signed zeros (never 0 * inf = NaN at a pole)
    const auto m = [s](T x) {
        return x == T() ? quat::detail::zero_prod(static_cast<T>(s), x) : static_cast<T>(static_cast<W>(x) * s);
    };
    return {static_cast<T>(r.real()), m(q.x()), m(q.y()), m(q.z())};
}
template <class T, class F>
[[nodiscard]] inline quaternion<T> lift(const quaternion<T>& q, F&& f, bool real_on_real) {
    using W = work_t<T>;
    return lift(q, f, real_on_real, [&f](W x) { return std::conj(f(cplx<W>(x, W()))); });
}

// ---- special values of the special functions (README, "Special functions") -------------
// A rule maps the i-axis image (x, y) of a special input to its value: y is the signed i
// component when the vector part is zero, |v| otherwise.  Rules are called for non-finite
// inputs without a NaN, and for finite inputs with a zero vector part (special points such
// as poles and signed zeros; any other point is left to the kernel).
namespace sv {
template <class T>
struct res {
    bool hit;
    std::complex<T> v;
};
template <class T>
[[nodiscard]] inline res<T> hit(T re, T im) { return {true, {re, im}}; }
template <class T>
[[nodiscard]] inline res<T> miss() { return {false, {}}; }
template <class T>
constexpr T inf = std::numeric_limits<T>::infinity();
template <class T>
[[nodiscard]] inline res<T> nn() { return hit(nan<T>(), nan<T>()); }
// limit of a component growing like c * infinity; a zero c is a zero limit (c's zero)
template <class T>
[[nodiscard]] inline T ci(T c) { return c == T() ? c : std::copysign(inf<T>, c); }
// the signed zero of y, +0 otherwise (not "y == 0 ? y : 0": MSVC folds that to +0)
template <class T>
[[nodiscard]] inline T zs(T y) { return y == T() ? std::copysign(T(), y) : T(); }
template <class T>
[[nodiscard]] inline res<T> cinf(T c, T s) { return hit(ci(c), ci(s)); }
// cos and sin of (p pi + e x), p with exact cospi/sinpi at half-integers, e = +-1
template <class T>
[[nodiscard]] inline std::complex<T> ph(T p, T e, T x) {
    const T cp = cospi(p), sp = sinpi(p);
    if (x == T()) return {cp, sp};
    const T cx = detail::lm::cos(x), sx = e * detail::lm::sin(x);
    return {cp * cx - sp * sx, sp * cx + cp * sx};
}
template <class T>
[[nodiscard]] inline res<T> cinf(const std::complex<T>& c) { return cinf(c.real(), c.imag()); }
// rays of a non-finite (x, y) with y >= +0: the real axis (y = +0), horizontal (y finite),
// vertical (x finite, y = +inf) and diagonal (both infinite)
enum kind { RP, RM, HP, HM, VV, DP, DM };
template <class T>
[[nodiscard]] inline kind ray(T x, T y) {
    if (std::isinf(y)) return std::isfinite(x) ? VV : x > T() ? DP : DM;
    if (y == T()) return x > T() ? RP : RM;
    return x > T() ? HP : HM;
}
// conjugate symmetry: a rule written for y >= +0 extended to y <= -0
template <class T, class R>
[[nodiscard]] inline res<T> csym(T x, T y, R&& r) {
    if (!std::signbit(y)) return r(x, y);
    res<T> o = r(x, -y);
    o.v = std::conj(o.v);
    return o;
}
template <class T>
[[nodiscard]] inline bool is_int(T x) { return std::isfinite(x) && x == std::floor(x); }
}  // namespace sv

// The lift with special values: finite inputs with a nonzero vector part take the plain lift
// (bit for bit), everything else goes through the rule; a NaN anywhere gives NaN (the zero
// vector part kept where f is real on the real axis).
template <class T, class F, class R, class G>
[[nodiscard]] inline quaternion<T> lift_s(const quaternion<T>& q, F&& f, bool real, R&& rule, G&& f_lower) {
    if (!quat::detail::needs_sv(q)) return lift(q, f, real, f_lower);
    if (quat::detail::is_finite(q)) {
        const sv::res<T> r = rule(q.w(), q.x());
        if (!r.hit) return lift(q, f, real, f_lower);
        return quat::detail::lift_sv(q, [&r](const std::complex<T>&) { return r.v; });
    }
    return quat::detail::lift_sv(q, [&](const std::complex<T>& z) -> std::complex<T> {
        const T a = z.real(), b = z.imag();
        if (std::isnan(a) || std::isnan(b)) return {nan<T>(), (b == T() && real) ? b : nan<T>()};
        const sv::res<T> r = rule(a, b);
        return r.hit ? r.v : std::complex<T>(nan<T>(), nan<T>());
    });
}
template <class T, class F, class R>
[[nodiscard]] inline quaternion<T> lift_s(const quaternion<T>& q, F&& f, bool real, R&& rule) {
    using W = work_t<T>;
    return lift_s(q, f, real, rule, [&f](W x) { return std::conj(f(cplx<W>(x, W()))); });
}

}  // namespace detail::sf

// ---- special-value rules, one per function (README table) ---------------------------------
// Each handles y >= +0 (conjugate symmetry gives y <= -0 where it holds).  "Z" in the
// comments is a zero limit whose sign is not specified (+0 is returned).
namespace detail::sf::sv {
#define QSV_FIN if (std::isfinite(x) && std::isfinite(y))
template <class T>
[[nodiscard]] inline res<T> r_tgamma(T x, T y) {
    return csym(x, y, [](T x, T y) -> res<T> {
        QSV_FIN {
            if (x == T()) return hit(std::copysign(inf<T>, x), y);
            if (x < T() && is_int(x)) return hit(nan<T>(), y);
            return miss<T>();
        }
        switch (ray(x, y)) {
            case RP: return hit(inf<T>, y);
            case RM: return hit(nan<T>(), y);
            case HP: case DP: return nn<T>();
            default: return hit(T(), T());  // |Gamma| -> 0
        }
    });
}
template <class T>
[[nodiscard]] inline res<T> r_lgamma(T x, T y) {
    return csym(x, y, [](T x, T y) -> res<T> {
        QSV_FIN return miss<T>();
        switch (ray(x, y)) {
            case RP: return hit(inf<T>, y);
            case RM: return hit(inf<T>, -inf<T>);  // C lgamma(-inf) = +inf; Im -> -inf above the cut
            case HP: case DP: return hit(inf<T>, inf<T>);
            case HM: return hit(-inf<T>, -inf<T>);
            default: return hit(-inf<T>, inf<T>);
        }
    });
}
template <class T>
[[nodiscard]] inline res<T> r_digamma(T x, T y) {
    return csym(x, y, [](T x, T y) -> res<T> {
        const T pi = quat::detail::pi_v<T>;
        QSV_FIN {
            if (x == T()) return hit(-std::copysign(inf<T>, x), y);
            if (x < T() && is_int(x)) return hit(nan<T>(), y);
            return miss<T>();
        }
        switch (ray(x, y)) {
            case RP: return hit(inf<T>, y);
            case RM: return hit(nan<T>(), y);
            case HP: return hit(inf<T>, T());
            case HM: return hit(inf<T>, nan<T>());
            case VV: return hit(inf<T>, pi / 2);
            case DP: return hit(inf<T>, pi / 4);
            default: return hit(inf<T>, 3 * pi / 4);
        }
    });
}
template <class T>
[[nodiscard]] inline res<T> r_polygamma(unsigned m, T x, T y) {
    if (m == 0) return r_digamma(x, y);
    return csym(x, y, [m](T x, T y) -> res<T> {
        QSV_FIN {
            if (x == T()) return hit(m % 2 ? inf<T> : -std::copysign(inf<T>, x), y);
            if (x < T() && is_int(x)) return hit(m % 2 ? inf<T> : nan<T>(), y);
            return miss<T>();
        }
        switch (ray(x, y)) {
            case RP: return hit(T(), y);
            case RM: return hit(nan<T>(), y);
            case HM: return nn<T>();
            default: return hit(T(), T());
        }
    });
}
template <class T>
[[nodiscard]] inline res<T> r_erf(T x, T y) {
    return csym(x, y, [](T x, T y) -> res<T> {
        QSV_FIN return x == T() ? hit(x, y) : miss<T>();
        switch (ray(x, y)) {
            case RP: return hit(T(1), y);
            case RM: return hit(T(-1), y);
            case HP: case DP: return hit(T(1), T());
            case HM: case DM: return hit(T(-1), T());
            default: return x == T() ? hit(x, inf<T>) : nn<T>();
        }
    });
}
template <class T>
[[nodiscard]] inline res<T> r_erfc(T x, T y) {
    return csym(x, y, [](T x, T y) -> res<T> {
        QSV_FIN return miss<T>();
        switch (ray(x, y)) {
            case RP: return hit(T(), y);
            case RM: return hit(T(2), y);
            case HP: case DP: return hit(T(), T());
            case HM: case DM: return hit(T(2), T());
            default: return x == T() ? hit(T(1), -inf<T>) : nn<T>();
        }
    });
}
template <class T>
[[nodiscard]] inline res<T> r_zeta(T x, T y) {
    return csym(x, y, [](T x, T y) -> res<T> {
        QSV_FIN return miss<T>();
        switch (ray(x, y)) {
            case RP: return hit(T(1), y);
            case RM: return hit(nan<T>(), y);
            case HP: case DP: return hit(T(1), T());
            default: return nn<T>();
        }
    });
}
template <class T>
[[nodiscard]] inline res<T> r_hurwitz(T a, T x, T y) {
    if (!(a > T())) return std::isfinite(x) && std::isfinite(y) ? miss<T>() : nn<T>();
    return csym(x, y, [a](T x, T y) -> res<T> {
        QSV_FIN return miss<T>();
        const kind k = ray(x, y);
        if (k == RM) return hit(nan<T>(), y);
        if (k != RP && k != HP && k != DP) return nn<T>();
        if (a == T(1)) return hit(T(1), k == RP ? y : T());
        if (a > T(1)) return hit(T(), k == RP ? y : T());
        // a^{-s} dominates: infinite with the phase -y log a
        if (k == RP) return hit(inf<T>, y);
        if (k == DP) return nn<T>();
        const T t = y * std::log(a);
        return cinf(detail::lm::cos(t), -detail::lm::sin(t));
    });
}
// Bessel functions: x finite and y -> +inf grows like e^y with the phase given
template <class T>
[[nodiscard]] inline T izero_sign(int n, T x) {  // sign of x^n for an integer power at a zero x
    return (std::signbit(x) && (n % 2 != 0)) ? T(-1) : T(1);
}
template <class T>
[[nodiscard]] inline res<T> r_bessel_j(T nu, T x, T y) {
    return csym(x, y, [nu](T x, T y) -> res<T> {
        QSV_FIN {
            if (x == T() && is_int(nu) && nu != T()) {
                const int n = static_cast<int>(std::abs(nu));  // J_{-n} = (-1)^n J_n
                const T s = izero_sign(n, x) * ((nu < T() && n % 2) ? T(-1) : T(1));
                return hit(std::copysign(T(), s), y);
            }
            return miss<T>();
        }
        switch (ray(x, y)) {
            case RP: case RM: return hit(T(), y);
            case HP: case HM: return hit(T(), T());
            case VV: return cinf(ph(nu / 2, T(-1), x));
            default: return nn<T>();
        }
    });
}
template <class T>
[[nodiscard]] inline res<T> r_bessel_y(T nu, T x, T y) {
    return csym(x, y, [nu](T x, T y) -> res<T> {
        QSV_FIN return miss<T>();
        switch (ray(x, y)) {
            case RP: case RM: return hit(T(), y);
            case HP: case HM: return hit(T(), T());
            case VV: return cinf(ph((nu + 1) / 2, T(-1), x));
            default: return nn<T>();
        }
    });
}
template <class T>
[[nodiscard]] inline res<T> r_bessel_i(T nu, T x, T y) {
    return csym(x, y, [nu](T x, T y) -> res<T> {
        QSV_FIN {
            if (x == T() && is_int(nu) && nu != T()) {  // I_{-n} = I_n
                const int n = static_cast<int>(std::abs(nu));
                return hit(std::copysign(T(), izero_sign(n, x)), y);
            }
            return miss<T>();
        }
        switch (ray(x, y)) {
            case RP: return hit(inf<T>, y);
            case RM: {
                const T c = cospi(nu), s = sinpi(nu);
                return hit(ci(c), s == T() ? y : ci(s));
            }
            case HP: return cinf(detail::lm::cos(y), detail::lm::sin(y));
            case HM: return cinf(ph(nu, T(-1), y));
            case VV: return hit(T(), T());
            default: return nn<T>();
        }
    });
}
template <class T>
[[nodiscard]] inline res<T> r_bessel_k(T, T x, T y) {
    return csym(x, y, [](T x, T y) -> res<T> {
        QSV_FIN return miss<T>();
        switch (ray(x, y)) {
            case RP: return hit(T(), y);
            case RM: return hit(T(), -inf<T>);
            case HM: return cinf(-detail::lm::sin(y), -detail::lm::cos(y));
            case DM: return nn<T>();
            default: return hit(T(), T());
        }
    });
}
// Hankel functions: H1 decays in the closed upper half plane, H2(z) = conj H1(conj z)
template <class T>
[[nodiscard]] inline res<T> r_hankel(int kind_, T nu, T x, T y) {
    QSV_FIN return miss<T>();
    if (kind_ == 1 || !std::isinf(y)) return hit(T(), T());
    return std::isfinite(x) ? cinf(ph(nu / 2, T(-1), x)) : nn<T>();
}
template <class T>
[[nodiscard]] inline res<T> r_sph_j(unsigned n, T x, T y) {
    return csym(x, y, [n](T x, T y) -> res<T> {
        QSV_FIN {
            if (x == T() && n > 0) return hit(std::copysign(T(), izero_sign(static_cast<int>(n % 2), x)), y);
            return miss<T>();
        }
        switch (ray(x, y)) {
            case RP: case RM: return hit(T(), y);
            case HP: case HM: return hit(T(), T());
            case VV: return cinf(ph(T(n) / 2, T(-1), x));
            default: return nn<T>();
        }
    });
}
template <class T>
[[nodiscard]] inline res<T> r_sph_y(unsigned n, T x, T y) {
    return csym(x, y, [n](T x, T y) -> res<T> {
        QSV_FIN {
            // y_n(x) ~ -(2n-1)!! / x^{n+1}
            if (x == T()) return hit(std::signbit(x) && n % 2 == 0 ? inf<T> : -inf<T>, y);
            return miss<T>();
        }
        switch (ray(x, y)) {
            case RP: case RM: return hit(T(), y);
            case HP: case HM: return hit(T(), T());
            case VV: return cinf(ph(T(n + 1) / 2, T(-1), x));
            default: return nn<T>();
        }
    });
}
template <class T>
[[nodiscard]] inline res<T> r_airy(bool bi, T x, T y) {
    return csym(x, y, [bi](T x, T y) -> res<T> {
        QSV_FIN return miss<T>();
        switch (ray(x, y)) {
            case RP: return hit(bi ? inf<T> : T(), y);
            case RM: return hit(T(), y);
            case HP: case DP: return bi ? nn<T>() : hit(T(), T());
            default: return nn<T>();
        }
    });
}
// elliptic integrals
template <class T>
[[nodiscard]] inline res<T> r_to_zero(T x, T y) {  // functions tending to 0 on every ray
    QSV_FIN return miss<T>();
    return hit(T(), zs(y));
}
template <class T>
[[nodiscard]] inline res<T> r_comp_ellint_2(T x, T y) {
    // E(k) = -i k + O(log k / k) for y >= +0
    return csym(x, y, [](T x, T y) -> res<T> {
        QSV_FIN return miss<T>();
        return hit(std::isinf(y) ? inf<T> : y, -x);
    });
}
template <class T>
[[nodiscard]] inline res<T> r_ellint_1(T k, T x, T y) {
    if (!(std::abs(k) <= T(1))) return std::isfinite(x) && std::isfinite(y) ? miss<T>() : nn<T>();
    return csym(x, y, [k](T x, T y) -> res<T> {
        QSV_FIN return miss<T>();
        if (k == T()) return hit(x, y);  // F(phi, 0) = phi
        using W = work_t<T>;
        switch (ray(x, y)) {
            case RP: case RM: return hit(x, y);
            case HP: case HM: return hit(x, nan<T>());  // Im is periodic in x
            default: {
                // F(x + i inf) = 2 K n + i K', n = round(x / pi), K' = K(sqrt(1 - k^2))
                const T kp = static_cast<T>(comp_ellint_1k(cplx<W>(static_cast<W>(detail::lm::sqrt((T(1) - k) * (T(1) + k))))).real());
                if (std::isinf(x)) return hit(x, kp);
                const T n = std::round(x / quat::detail::pi_v<T>);
                if (n == T()) return hit(T(), kp);
                const T K = static_cast<T>(comp_ellint_1k(cplx<W>(static_cast<W>(k))).real());
                return hit(2 * K * n, kp);
            }
        }
    });
}
template <class T>
[[nodiscard]] inline res<T> r_ellint_2(T k, T x, T y) {
    if (!(std::abs(k) <= T(1))) return std::isfinite(x) && std::isfinite(y) ? miss<T>() : nn<T>();
    return csym(x, y, [k](T x, T y) -> res<T> {
        QSV_FIN return miss<T>();
        if (k == T()) return hit(x, y);  // E(phi, 0) = phi
        switch (ray(x, y)) {
            case RP: case RM: return hit(x, y);
            case HP: case HM: return hit(x, nan<T>());
            case VV: {
                // E ~ (-1)^n (i k / 2) e^{-i (x + i y)}, n = round(x / pi)
                const T n = std::round(x / quat::detail::pi_v<T>);
                const T s = std::fmod(n, T(2)) == T() ? T(1) : T(-1);
                return cinf(s * detail::lm::sin(x), s * detail::lm::cos(x));
            }
            default: return nn<T>();
        }
    });
}
// exponential, sine/cosine and Fresnel integrals
template <class T>
[[nodiscard]] inline res<T> r_ei(T x, T y) {
    return csym(x, y, [](T x, T y) -> res<T> {
        const T pi = quat::detail::pi_v<T>;
        QSV_FIN return miss<T>();
        switch (ray(x, y)) {
            case RP: return hit(inf<T>, y);
            case RM: return hit(-T(), y);
            case HP: return cinf(detail::lm::cos(y), detail::lm::sin(y));
            case DP: return nn<T>();
            default: return hit(T(), pi);
        }
    });
}
template <class T>
[[nodiscard]] inline res<T> r_en(unsigned n, T x, T y) {
    return csym(x, y, [n](T x, T y) -> res<T> {
        QSV_FIN return miss<T>();
        switch (ray(x, y)) {
            case RP: return hit(T(), y);
            case RM: return n == 0 ? hit(-inf<T>, y) : hit(-inf<T>, n == 1 ? -quat::detail::pi_v<T> : -inf<T>);
            case HM: return cinf(-detail::lm::cos(y), detail::lm::sin(y));
            case DM: return nn<T>();
            default: return hit(T(), T());
        }
    });
}
template <class T>
[[nodiscard]] inline res<T> r_si(T x, T y) {
    return csym(x, y, [](T x, T y) -> res<T> {
        const T h = quat::detail::pi_v<T> / 2;
        QSV_FIN return x == T() ? hit(x, y) : miss<T>();
        switch (ray(x, y)) {
            case RP: return hit(h, y);
            case RM: return hit(-h, y);
            case HP: return hit(h, T());
            case HM: return hit(-h, T());
            case VV: return cinf(detail::lm::sin(x), detail::lm::cos(x));
            default: return nn<T>();
        }
    });
}
template <class T>
[[nodiscard]] inline res<T> r_ci(T x, T y) {
    return csym(x, y, [](T x, T y) -> res<T> {
        const T pi = quat::detail::pi_v<T>;
        QSV_FIN return x == T() ? hit(-inf<T>, std::signbit(x) ? pi : y) : miss<T>();
        switch (ray(x, y)) {
            case RP: return hit(T(), y);
            case HP: return hit(T(), T());
            case RM: case HM: return hit(T(), pi);
            case VV: return x == T() ? hit(inf<T>, pi / 2) : cinf(detail::lm::cos(x), -detail::lm::sin(x));
            default: return nn<T>();
        }
    });
}
template <class T>
[[nodiscard]] inline res<T> r_fresnel(bool s, T x, T y) {
    return csym(x, y, [s](T x, T y) -> res<T> {
        QSV_FIN return x == T() ? hit(x, y) : miss<T>();
        switch (ray(x, y)) {
            case RP: return hit(T(0.5), y);
            case RM: return hit(T(-0.5), y);
            case VV: return x == T() ? hit(x, s ? T(-0.5) : T(0.5)) : nn<T>();
            case DP: return s ? hit(-inf<T>, inf<T>) : hit(inf<T>, inf<T>);
            case DM: return s ? hit(inf<T>, inf<T>) : hit(-inf<T>, inf<T>);
            default: return nn<T>();
        }
    });
}
// polylogarithm: Li_s(z) ~ -log(-z)^s / Gamma(s + 1)
template <class T>
[[nodiscard]] inline res<T> r_polylog(T s, T x, T y) {
    return csym(x, y, [s](T x, T y) -> res<T> {
        const T pi = quat::detail::pi_v<T>;
        QSV_FIN return x == T() ? hit(x, y) : miss<T>();
        const kind k = ray(x, y);
        const bool real = k == RM;
        if (!(s > T())) return hit(s == T() ? T(-1) : T(), real ? y : T());
        if (k == RM || k == HM) return hit(-inf<T>, real ? y : T());
        const T th = (k == RP || k == HP) ? -pi : k == VV ? -pi / 2 : k == DP ? -3 * pi / 4 : -pi / 4;
        return hit(-inf<T>, s > T(1) ? inf<T> : s == T(1) ? -th : T());
    });
}
// Lambert W_k: W -> +inf + i (arg z + 2 pi k) on every ray
template <class T>
[[nodiscard]] inline res<T> r_lambert_w(int kb, T x, T y) {
    const T pi = quat::detail::pi_v<T>;
    QSV_FIN return x == T() && kb == 0 ? hit(x, y) : miss<T>();
    T th;
    if (std::isinf(x) && std::isinf(y)) th = std::copysign(x > T() ? pi / 4 : 3 * pi / 4, y);
    else if (std::isinf(y)) th = std::copysign(pi / 2, y);
    else if (x > T()) th = T();
    else th = std::copysign(pi, y);
    const T im = th + 2 * pi * T(kb);
    return hit(inf<T>, im == T() ? zs(y) : im);
}
// incomplete gamma: Gamma(a, z) ~ z^{a-1} e^{-z}
template <class T>
[[nodiscard]] inline res<T> r_gamma_inc(T a, bool upper, bool reg, T x, T y) {
    if (reg && !(a > T())) return std::isfinite(x) && std::isfinite(y) ? miss<T>() : nn<T>();
    if (!upper && a <= T() && is_int(a)) return std::isfinite(x) && std::isfinite(y) ? miss<T>() : nn<T>();
    return csym(x, y, [a, upper, reg](T x, T y) -> res<T> {
        QSV_FIN return miss<T>();
        // the constant of the lower function: gamma(a, inf) = Gamma(a), P(a, inf) = 1
        const T g = reg ? T(1) : static_cast<T>(detail::lm::tgamma(a));
        const T top = upper ? T() : g;
        const T sg = upper ? T(1) : T(-1);  // gamma = Gamma(a) - Gamma(a, z)
        switch (ray(x, y)) {
            case RP: return hit(top, y);
            case HP: case DP: return hit(top, T());
            case VV: return a < T(1) ? hit(top, T()) : nn<T>();
            case DM: return nn<T>();
            case HM: {
                const std::complex<T> c = ph(a - 1, T(-1), y);
                return cinf(sg * c.real(), sg * c.imag());
            }
            default: {  // RM, above the cut
                if (is_int(a) && a > T()) return hit(ci(sg * cospi(a - 1)), y);
                if (is_int(a)) {  // Gamma(-n, z): (-1)^n / n! (E1(z) - ...), Im -> (-1)^{n+1} pi / n!
                    const T n = -a, s = std::fmod(n, T(2)) == T() ? T(-1) : T(1);
                    return hit(s * inf<T>, s * quat::detail::pi_v<T> / static_cast<T>(detail::lm::tgamma(n + 1)));
                }
                const T c = cospi(a), s = sinpi(a);
                const T re = c == T() ? (upper ? (reg ? T(1) : g) : T()) : ci(-sg * c);
                return hit(re, ci(-sg * s));
            }
        }
    });
}
template <class T>
[[nodiscard]] inline res<T> r_beta(T b, T x, T y) {
    if (!(b > T())) return std::isfinite(x) && std::isfinite(y) ? miss<T>() : nn<T>();
    return csym(x, y, [b](T x, T y) -> res<T> {
        QSV_FIN return miss<T>();
        const kind k = ray(x, y);
        if (k == RM) return hit(is_int(b) ? T() : nan<T>(), y);
        return hit(T(), k == RP ? y : T());
    });
}

// Real polynomials sum c_k z^k: the limit along the ray of a non-finite (x, y), y >= +0.
template <class T>
[[nodiscard]] inline std::complex<T> poly_ray(const std::vector<long double>& c, T x, T y) {
    using L = long double;
    int n = static_cast<int>(c.size()) - 1;
    while (n > 0 && c[static_cast<std::size_t>(n)] == L(0)) --n;
    const auto C = [&c](int k) { return c[static_cast<std::size_t>(k)]; };
    const auto lim = [](int k, L v) { return k == 0 ? static_cast<T>(v) : std::copysign(inf<T>, static_cast<T>(v)); };
    if (n == 0) return {static_cast<T>(C(0)), zs(y)};
    if (std::isfinite(y)) {  // x = +-inf: Re ~ c_n x^n, Im ~ n c_n x^{n-1} y
        const L sx = x > T() ? L(1) : L(-1);
        const L re = C(n) * (n % 2 ? sx : L(1));
        if (y == T()) return {lim(1, re), y};
        if (n == 1) return {lim(1, re), static_cast<T>(C(1) * static_cast<L>(y))};
        return {lim(1, re), lim(1, C(n) * ((n - 1) % 2 ? sx : L(1)) * static_cast<L>(y))};
    }
    if (std::isfinite(x)) {  // y = +inf: p(x + i t) = sum a_k (i t)^k, a_k the Taylor coefficients at x
        std::vector<L> a(c.begin(), c.begin() + n + 1);
        const L lx = static_cast<L>(x);
        for (int i = 0; i < n; ++i)
            for (int j = n - 1; j >= i; --j) a[static_cast<std::size_t>(j)] += lx * a[static_cast<std::size_t>(j + 1)];
        T re = T(), im = T();
        for (int k = n; k >= 0; --k)
            if (k % 2 == 0 && a[static_cast<std::size_t>(k)] != L(0)) { re = lim(k, a[static_cast<std::size_t>(k)] * ((k / 2) % 2 ? L(-1) : L(1))); break; }
        for (int k = n; k >= 1; --k)
            if (k % 2 == 1 && a[static_cast<std::size_t>(k)] != L(0)) { im = lim(k, a[static_cast<std::size_t>(k)] * (((k - 1) / 2) % 2 ? L(-1) : L(1))); break; }
        return {re, im};
    }
    // diagonal z = t (sx + i): coefficients c_k (sx + i)^k, exact Gaussian integers
    const L sx = x > T() ? L(1) : L(-1);
    std::vector<L> dr(static_cast<std::size_t>(n) + 1), di(static_cast<std::size_t>(n) + 1);
    dr[0] = 1; di[0] = 0;
    for (int k = 1; k <= n; ++k) {
        const L r = dr[static_cast<std::size_t>(k - 1)], m = di[static_cast<std::size_t>(k - 1)];
        dr[static_cast<std::size_t>(k)] = r * sx - m;
        di[static_cast<std::size_t>(k)] = r + m * sx;
    }
    T re = T(), im = T();
    for (int k = n; k >= 0; --k)
        if (C(k) * dr[static_cast<std::size_t>(k)] != L(0)) { re = lim(k, C(k) * dr[static_cast<std::size_t>(k)]); break; }
    for (int k = n; k >= 1; --k)
        if (C(k) * di[static_cast<std::size_t>(k)] != L(0)) { im = lim(k, C(k) * di[static_cast<std::size_t>(k)]); break; }
    return {re, im};
}
// monomial coefficients of the orthogonal polynomials, from their three-term recurrences
// p_{k+1} = (A_k x + B_k) p_k - C_k p_{k-1}
template <class F>
[[nodiscard]] inline std::vector<long double> poly_rec(unsigned n, std::vector<long double> p0, std::vector<long double> p1, F&& abc) {
    using L = long double;
    if (n == 0) return p0;
    for (unsigned k = 1; k < n; ++k) {
        L A, B, Cc;
        abc(k, A, B, Cc);
        std::vector<L> p2(p1.size() + 1, L(0));
        for (std::size_t j = 0; j < p1.size(); ++j) {
            p2[j + 1] += A * p1[j];
            p2[j] += B * p1[j];
        }
        for (std::size_t j = 0; j < p0.size(); ++j) p2[j] -= Cc * p0[j];
        p0 = std::move(p1);
        p1 = std::move(p2);
    }
    return p1;
}
[[nodiscard]] inline std::vector<long double> poly_deriv(std::vector<long double> c, unsigned m) {
    for (unsigned i = 0; i < m && !c.empty(); ++i) {
        for (std::size_t j = 1; j < c.size(); ++j) c[j - 1] = c[j] * static_cast<long double>(j);
        c.pop_back();
    }
    if (c.empty()) c.push_back(0);
    return c;
}
[[nodiscard]] inline std::vector<long double> poly_mul(const std::vector<long double>& a, const std::vector<long double>& b) {
    std::vector<long double> r(a.size() + b.size() - 1, 0);
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < b.size(); ++j) r[i + j] += a[i] * b[j];
    return r;
}
// Special values of a polynomial-type function evaluated in quaternion arithmetic: rule(x, y)
// gives the value for a non-finite (x, y) (y may be a signed zero); NaN as in lift_s.
template <class T, class R>
[[nodiscard]] inline quaternion<T> poly_sv(const quaternion<T>& q, bool real, R&& rule) {
    return quat::detail::lift_sv(q, [&](const std::complex<T>& z) -> std::complex<T> {
        const T a = z.real(), b = z.imag();
        if (std::isnan(a) || std::isnan(b)) return {nan<T>(), (b == T() && real) ? b : nan<T>()};
        return rule(a, b);
    });
}
template <class T>
[[nodiscard]] inline std::complex<T> poly_ray_sym(const std::vector<long double>& c, T x, T y) {
    if (!std::signbit(y)) return poly_ray(c, x, y);
    return std::conj(poly_ray(c, x, -y));
}
#undef QSV_FIN
}  // namespace detail::sf::sv

#define QUAT_SF_LIFT(NAME, EXPR, REAL)                                                   \
    template <class T>                                                                   \
    [[nodiscard]] inline quaternion<T> NAME(const quaternion<T>& q) noexcept {           \
        return detail::sf::lift(                                                         \
            q, [](const auto& z) { return EXPR; }, REAL);                                \
    }

// Inverse trigonometric and hyperbolic functions (std::complex branch cuts).  Non-finite
// inputs and inputs with a zero vector part (signed zero on the cut) follow Annex G
// through quat::detail::lift_sv and std::complex.
namespace detail::sf::ag {
// Annex G (C11 G.6.1-G.6.2) kernels for the inverse functions on std::complex<T>, reduced to
// y >= +0 (conjugate symmetry) and, for the odd functions, x >= +0; finite points use the
// accurate kernels above.
using quat::detail::ag::negy;
using quat::detail::ag::cj;
using quat::detail::ag::nn;
using quat::detail::ag::mi;
template <class T>
using C = std::complex<T>;
// finite kernel at x + i y with zeros of the result given fixed signs (zre, zim)
template <class T, class F>
[[nodiscard]] inline C<T> fk(F f, T x, T y, T zre, T zim) {
    using W = work_t<T>;
    const cplx<W> r = f(cplx<W>(static_cast<W>(x), static_cast<W>(y)));
    const T re = static_cast<T>(r.real()), im = static_cast<T>(r.imag());
    return {re == T() ? zre : re, im == T() ? zim : im};
}
template <class T>
[[nodiscard]] inline C<T> asinh(T x, T y) {
    if (!std::isnan(x) && std::signbit(x)) return -asinh(-x, -y);
    if (negy(y)) return cj(asinh(x, -y));
    const T inf = std::numeric_limits<T>::infinity(), pi = quat::detail::pi_v<T>;
    if (std::isinf(y)) return std::isnan(x) ? C<T>(inf, x) : C<T>(inf, std::isinf(x) ? pi / 4 : pi / 2);
    if (std::isinf(x)) return {inf, std::isnan(y) ? y : T()};
    if (std::isnan(x)) return y == T() ? C<T>(x, y) : nn<T>();
    if (std::isnan(y)) return nn<T>();
    if (x == T() && y == T()) return {x, y};
    return fk([](const auto& z) { return sf::asinh(z); }, x, y, T(), T());
}
template <class T>
[[nodiscard]] inline C<T> atanh(T x, T y) {
    if (!std::isnan(x) && std::signbit(x)) return -atanh(-x, -y);
    if (negy(y)) return cj(atanh(x, -y));
    const T inf = std::numeric_limits<T>::infinity(), pi = quat::detail::pi_v<T>;
    if (std::isinf(y)) return {T(), pi / 2};
    if (std::isinf(x)) return {T(), std::isnan(y) ? y : pi / 2};
    if (std::isnan(x)) return nn<T>();
    if (std::isnan(y)) return x == T() ? C<T>(x, y) : nn<T>();
    if (y == T()) {
        if (x == T()) return {x, y};
        if (x == T(1)) return {inf, y};
    }
    return fk([](const auto& z) { return sf::atanh(z); }, x, y, T(), T());
}
template <class T>
[[nodiscard]] inline C<T> acosh(T x, T y) {
    if (negy(y)) return cj(acosh(x, -y));
    const T inf = std::numeric_limits<T>::infinity(), pi = quat::detail::pi_v<T>;
    if (std::isinf(y))
        return std::isnan(x) ? C<T>(inf, x)
                             : C<T>(inf, x == -inf ? 3 * pi / 4 : x == inf ? pi / 4 : pi / 2);
    if (std::isnan(y)) return std::isinf(x) ? C<T>(inf, y) : nn<T>();
    if (std::isinf(x)) return {inf, x < T() ? pi : T()};
    if (std::isnan(x)) return nn<T>();
    if (x == T() && y == T()) return {T(), pi / 2};
    return fk([](const auto& z) { return sf::acosh(z); }, x, y, T(), T());
}
template <class T>
[[nodiscard]] inline C<T> acos(T x, T y) {
    if (negy(y)) return cj(acos(x, -y));
    const T inf = std::numeric_limits<T>::infinity(), pi = quat::detail::pi_v<T>;
    if (std::isinf(y))
        return std::isnan(x) ? C<T>(x, -inf)
                             : C<T>(x == -inf ? 3 * pi / 4 : x == inf ? pi / 4 : pi / 2, -inf);
    if (std::isnan(y)) return std::isinf(x) ? C<T>(y, inf) : x == T() ? C<T>(pi / 2, y) : nn<T>();
    if (std::isinf(x)) return {x < T() ? pi : T(), -inf};
    if (std::isnan(x)) return nn<T>();
    if (x == T() && y == T()) return {pi / 2, -T()};
    return fk([](const auto& z) { return sf::acos(z); }, x, y, T(), -T());
}
template <class T>
[[nodiscard]] inline C<T> asin(T x, T y) { return mi(asinh(-y, x)); }
template <class T>
[[nodiscard]] inline C<T> atan(T x, T y) { return mi(atanh(-y, x)); }
}  // namespace detail::sf::ag
#define QUAT_SF_LIFT_G(NAME, EXPR)                                                       \
    template <class T>                                                                   \
    [[nodiscard]] inline quaternion<T> NAME(const quaternion<T>& q) noexcept {           \
        if (quat::detail::needs_sv(q))                                                   \
            return quat::detail::lift_sv(                                                \
                q, [](const std::complex<T>& z) { return detail::sf::ag::NAME(z.real(), z.imag()); }); \
        return detail::sf::lift(                                                         \
            q, [](const auto& z) { return EXPR; }, false);                               \
    }
QUAT_SF_LIFT_G(asin, detail::sf::asin(z))
QUAT_SF_LIFT_G(acos, detail::sf::acos(z))
QUAT_SF_LIFT_G(atan, detail::sf::atan(z))
QUAT_SF_LIFT_G(asinh, detail::sf::asinh(z))
QUAT_SF_LIFT_G(acosh, detail::sf::acosh(z))
QUAT_SF_LIFT_G(atanh, detail::sf::atanh(z))
#undef QUAT_SF_LIFT_G

// Elementary functions with good accuracy near zero.
QUAT_SF_LIFT(expm1, detail::sf::expm1(z), true)
QUAT_SF_LIFT(log1p, detail::sf::log1p(z), false)
QUAT_SF_LIFT(cbrt, detail::sf::cbrt(z), false)

#undef QUAT_SF_LIFT

// Gamma family, error functions, Riemann zeta (special values: detail::sf::sv rules).
#define QUAT_SF_LIFT_S(NAME, EXPR, REAL, RULE)                                           \
    template <class T>                                                                   \
    [[nodiscard]] inline quaternion<T> NAME(const quaternion<T>& q) noexcept {           \
        return detail::sf::lift_s(                                                       \
            q, [](const auto& z) { return EXPR; }, REAL,                                 \
            [](T x, T y) { return detail::sf::sv::RULE(x, y); });                        \
    }
QUAT_SF_LIFT_S(tgamma, detail::sf::tgamma(z), true, r_tgamma)
QUAT_SF_LIFT_S(lgamma, detail::sf::lgamma(z), false, r_lgamma)
QUAT_SF_LIFT_S(digamma, detail::sf::digamma(z), true, r_digamma)
QUAT_SF_LIFT_S(erf, detail::sf::erf(z), true, r_erf)
QUAT_SF_LIFT_S(erfc, detail::sf::erfc(z), true, r_erfc)
QUAT_SF_LIFT_S(riemann_zeta, detail::sf::riemann_zeta(z), true, r_zeta)
#undef QUAT_SF_LIFT_S

/// Cylinder Bessel function of the first kind J_n(q) for integer order n (entire in q).
template <class T>
[[nodiscard]] inline quaternion<T> cyl_bessel_j(int n, const quaternion<T>& q) noexcept {
    return detail::sf::lift_s(q, [n](const auto& z) { return detail::sf::cyl_bessel_j(n, z); }, true,
                              [n](T x, T y) { return detail::sf::sv::r_bessel_j(T(n), x, y); });
}

/// Euler beta function B(q, b) = Gamma(q) Gamma(b) / Gamma(q + b) for real b.
/// (For two quaternions on different axes the product is not well defined.)
template <class T>
[[nodiscard]] inline quaternion<T> beta(const quaternion<T>& q, const T& b) noexcept {
    using W = detail::sf::work_t<T>;
    const W wb = static_cast<W>(b);
    return detail::sf::lift_s(
        q,
        [wb](const std::complex<W>& z) {
            using namespace detail::sf;
            const cplx<W> cb(wb);
            if (is_nonpositive_integer(z) || is_nonpositive_integer(cb))
                return cplx<W>(nan<W>(), W());
            return detail::lm::exp(lgamma(z) + lgamma(cb) - lgamma(z + cb));
        },
        true, [b](T x, T y) { return detail::sf::sv::r_beta(b, x, y); });
}
template <class T>
[[nodiscard]] inline quaternion<T> beta(const T& a, const quaternion<T>& q) noexcept {
    return beta(q, a);
}


// =====================================================================================
// Extended special functions (public API)
// =====================================================================================
// Real parameters (orders, degrees, moduli) are taken as a separate floating type N or an
// integer; the quaternion argument is lifted by axis reduction.  The `pred` passed to the
// lift tells, for a purely real argument x, whether the function is real there.

namespace detail::sf {
template <class T, class F, class P>
    requires std::is_invocable_r_v<bool, P&, work_t<T>>
[[nodiscard]] inline quaternion<T> lift(const quaternion<T>& q, F&& f, P&& real_at) {
    using W = work_t<T>;
    if (q.x() == T() && q.y() == T() && q.z() == T()) return lift(q, f, static_cast<bool>(real_at(static_cast<W>(q.w()))));
    return lift(q, f, false);
}
template <class T, class F, class P, class R, class G>
    requires std::is_invocable_r_v<bool, P&, work_t<T>>
[[nodiscard]] inline quaternion<T> lift_s(const quaternion<T>& q, F&& f, P&& real_at, R&& rule, G&& f_lower) {
    using W = work_t<T>;
    const bool real = q.x() == T() && q.y() == T() && q.z() == T() && real_at(static_cast<W>(q.w()));
    return lift_s(q, f, real, rule, f_lower);
}
template <class T, class F, class P, class R>
    requires std::is_invocable_r_v<bool, P&, work_t<T>>
[[nodiscard]] inline quaternion<T> lift_s(const quaternion<T>& q, F&& f, P&& real_at, R&& rule) {
    using W = work_t<T>;
    const bool real = q.x() == T() && q.y() == T() && q.z() == T() && real_at(static_cast<W>(q.w()));
    return lift_s(q, f, real, rule);
}
template <class N, class T>
[[nodiscard]] inline work_t<T> par(N v) noexcept { return static_cast<work_t<T>>(v); }
}  // namespace detail::sf

#define QUAT_SF_W using W = detail::sf::work_t<T>
#define QUAT_SF_C std::complex<W>

// ---- Bessel functions of real order ------------------------------------------------------
/// J_nu(q), real order nu; principal branch with the cut along the negative real axis.
template <class T, std::floating_point N>
[[nodiscard]] inline quaternion<T> cyl_bessel_j(N nu, const quaternion<T>& q) {
    QUAT_SF_W;
    const W m = static_cast<W>(nu);
    return detail::sf::lift_s(q, [m](const QUAT_SF_C& z) { return detail::sf::bessel_j_nu(m, z); },
                            [m](W x) { return x >= W() || detail::sf::is_integer(m); },
                              [nu](T x, T y) { return detail::sf::sv::r_bessel_j(static_cast<T>(nu), x, y); });
}
/// Y_nu(q) (Bessel function of the second kind / Neumann function).
template <class T, std::floating_point N>
[[nodiscard]] inline quaternion<T> cyl_neumann(N nu, const quaternion<T>& q) {
    QUAT_SF_W;
    const W m = static_cast<W>(nu);
    return detail::sf::lift_s(q, [m](const QUAT_SF_C& z) { return detail::sf::bessel_y_nu(m, z); },
                            [](W x) { return x > W(); },
                              [nu](T x, T y) { return detail::sf::sv::r_bessel_y(static_cast<T>(nu), x, y); });
}
template <class T>
[[nodiscard]] inline quaternion<T> cyl_neumann(int n, const quaternion<T>& q) {
    return cyl_neumann(static_cast<double>(n), q);
}
/// I_nu(q) (modified Bessel function of the first kind).
template <class T, std::floating_point N>
[[nodiscard]] inline quaternion<T> cyl_bessel_i(N nu, const quaternion<T>& q) {
    QUAT_SF_W;
    const W m = static_cast<W>(nu);
    return detail::sf::lift_s(q, [m](const QUAT_SF_C& z) { return detail::sf::bessel_i_nu(m, z); },
                            [m](W x) { return x >= W() || detail::sf::is_integer(m); },
                              [nu](T x, T y) { return detail::sf::sv::r_bessel_i(static_cast<T>(nu), x, y); });
}
template <class T>
[[nodiscard]] inline quaternion<T> cyl_bessel_i(int n, const quaternion<T>& q) {
    return cyl_bessel_i(static_cast<double>(n), q);
}
/// K_nu(q) (modified Bessel function of the second kind).
template <class T, std::floating_point N>
[[nodiscard]] inline quaternion<T> cyl_bessel_k(N nu, const quaternion<T>& q) {
    QUAT_SF_W;
    const W m = static_cast<W>(nu);
    return detail::sf::lift_s(q, [m](const QUAT_SF_C& z) { return detail::sf::bessel_k_nu(m, z); },
                            [](W x) { return x > W(); },
                              [nu](T x, T y) { return detail::sf::sv::r_bessel_k(static_cast<T>(nu), x, y); });
}
template <class T>
[[nodiscard]] inline quaternion<T> cyl_bessel_k(int n, const quaternion<T>& q) {
    return cyl_bessel_k(static_cast<double>(n), q);
}
/// Hankel functions H1_nu(q) = J + iY and H2_nu(q) = J - iY ("i" is the axis of q).
template <class T, std::floating_point N>
[[nodiscard]] inline quaternion<T> cyl_hankel_1(N nu, const quaternion<T>& q) {
    QUAT_SF_W;
    const W m = static_cast<W>(nu);
    return detail::sf::lift_s(
        q, [m](const QUAT_SF_C& z) { return detail::sf::hankel_1(m, z); }, false,
        [nu](T x, T y) { return detail::sf::sv::r_hankel(1, static_cast<T>(nu), x, y); },
        [m](W x) { return std::conj(detail::sf::hankel_2(m, QUAT_SF_C(x, W()))); });
}
template <class T>
[[nodiscard]] inline quaternion<T> cyl_hankel_1(int n, const quaternion<T>& q) {
    return cyl_hankel_1(static_cast<double>(n), q);
}
template <class T, std::floating_point N>
[[nodiscard]] inline quaternion<T> cyl_hankel_2(N nu, const quaternion<T>& q) {
    QUAT_SF_W;
    const W m = static_cast<W>(nu);
    return detail::sf::lift_s(
        q, [m](const QUAT_SF_C& z) { return detail::sf::hankel_2(m, z); }, false,
        [nu](T x, T y) { return detail::sf::sv::r_hankel(2, static_cast<T>(nu), x, y); },
        [m](W x) { return std::conj(detail::sf::hankel_1(m, QUAT_SF_C(x, W()))); });
}
template <class T>
[[nodiscard]] inline quaternion<T> cyl_hankel_2(int n, const quaternion<T>& q) {
    return cyl_hankel_2(static_cast<double>(n), q);
}
/// Spherical Bessel functions j_n, y_n and Hankel functions h1_n, h2_n.
template <class T>
[[nodiscard]] inline quaternion<T> sph_bessel(unsigned n, const quaternion<T>& q) {
    QUAT_SF_W;
    return detail::sf::lift_s(q, [n](const QUAT_SF_C& z) { return detail::sf::sph_bessel_j(n, z); }, true,
                              [n](T x, T y) { return detail::sf::sv::r_sph_j(n, x, y); });
}
template <class T>
[[nodiscard]] inline quaternion<T> sph_neumann(unsigned n, const quaternion<T>& q) {
    QUAT_SF_W;
    return detail::sf::lift_s(q, [n](const QUAT_SF_C& z) { return detail::sf::sph_bessel_y(n, z); }, true,
                              [n](T x, T y) { return detail::sf::sv::r_sph_y(n, x, y); });
}
template <class T>
[[nodiscard]] inline quaternion<T> sph_hankel_1(unsigned n, const quaternion<T>& q) {
    QUAT_SF_W;
    return detail::sf::lift_s(
        q, [n](const QUAT_SF_C& z) { return detail::sf::sph_hankel(n, z, 1); }, false,
        [n](T x, T y) { return detail::sf::sv::r_hankel(1, static_cast<T>(n), x, y); },
        [n](W x) { return std::conj(detail::sf::sph_hankel(n, QUAT_SF_C(x, W()), 2)); });
}
template <class T>
[[nodiscard]] inline quaternion<T> sph_hankel_2(unsigned n, const quaternion<T>& q) {
    QUAT_SF_W;
    return detail::sf::lift_s(
        q, [n](const QUAT_SF_C& z) { return detail::sf::sph_hankel(n, z, 2); }, false,
        [n](T x, T y) { return detail::sf::sv::r_hankel(2, static_cast<T>(n), x, y); },
        [n](W x) { return std::conj(detail::sf::sph_hankel(n, QUAT_SF_C(x, W()), 1)); });
}
/// Airy functions Ai(q), Bi(q) (entire).
template <class T>
[[nodiscard]] inline quaternion<T> airy_ai(const quaternion<T>& q) {
    QUAT_SF_W;
    return detail::sf::lift_s(q, [](const QUAT_SF_C& z) { return detail::sf::airy(z, false); }, true,
                              [](T x, T y) { return detail::sf::sv::r_airy(false, x, y); });
}
template <class T>
[[nodiscard]] inline quaternion<T> airy_bi(const quaternion<T>& q) {
    QUAT_SF_W;
    return detail::sf::lift_s(q, [](const QUAT_SF_C& z) { return detail::sf::airy(z, true); }, true,
                              [](T x, T y) { return detail::sf::sv::r_airy(true, x, y); });
}

// ---- orthogonal polynomials (evaluated directly in quaternion arithmetic) -----------------
// All are polynomials in q with real coefficients, hence exact axis reductions.
namespace detail::sf {
template <class T>
using qw = quaternion<work_t<T>>;
template <class T>
[[nodiscard]] inline qw<T> to_w(const quaternion<T>& q) {
    using W = work_t<T>;
    return {static_cast<W>(q.w()), static_cast<W>(q.x()), static_cast<W>(q.y()), static_cast<W>(q.z())};
}
template <class T, class W>
[[nodiscard]] inline quaternion<T> from_w(const quaternion<W>& q) {
    return {static_cast<T>(q.w()), static_cast<T>(q.x()), static_cast<T>(q.y()), static_cast<T>(q.z())};
}
}  // namespace detail::sf

namespace detail::sf::sv {
using pc = std::vector<long double>;
[[nodiscard]] inline pc c_legendre(unsigned l) {
    return poly_rec(l, {1}, {0, 1}, [](unsigned k, long double& A, long double& B, long double& C) {
        A = (2.0L * k + 1) / (k + 1); B = 0; C = static_cast<long double>(k) / (k + 1);
    });
}
[[nodiscard]] inline pc c_hermite(unsigned n) {
    return poly_rec(n, {1}, {0, 2}, [](unsigned k, long double& A, long double& B, long double& C) {
        A = 2; B = 0; C = 2.0L * k;
    });
}
[[nodiscard]] inline pc c_laguerre(unsigned n, unsigned m) {
    const long double a = m;
    return poly_rec(n, {1}, {1 + a, -1}, [a](unsigned k, long double& A, long double& B, long double& C) {
        A = -1.0L / (k + 1); B = (2.0L * k + 1 + a) / (k + 1); C = (k + a) / (k + 1);
    });
}
[[nodiscard]] inline pc c_chebyshev(unsigned n, bool second) {
    return poly_rec(n, {1}, second ? pc{0, 2} : pc{0, 1}, [](unsigned, long double& A, long double& B, long double& C) {
        A = 2; B = 0; C = 1;
    });
}
[[nodiscard]] inline pc c_one_minus_x2(unsigned p) {  // (1 - x^2)^p
    pc r{1};
    for (unsigned i = 0; i < p; ++i) r = poly_mul(r, pc{1, 0, -1});
    return r;
}
// I_z(a, b) for positive integers a, b: sum_{j=a}^{n} C(n, j) z^j (1 - z)^{n-j}, n = a + b - 1
[[nodiscard]] inline pc c_ibeta(unsigned a, unsigned b) {
    const unsigned n = a + b - 1;
    pc r(n + 1, 0);
    long double binom = 1;  // C(n, j)
    for (unsigned j = 0; j <= n; ++j) {
        if (j >= a) {
            pc t{binom};
            for (unsigned i = 0; i < j; ++i) t = poly_mul(t, pc{0, 1});
            for (unsigned i = 0; i < n - j; ++i) t = poly_mul(t, pc{1, -1});
            for (std::size_t i = 0; i < t.size(); ++i) r[i] += t[i];
        }
        binom = binom * (n - j) / (j + 1);
    }
    return r;
}
template <class T>
[[nodiscard]] inline res<T> r_ibeta(T a, T b, bool comp, T x, T y) {
    if (std::isfinite(x) && std::isfinite(y)) return miss<T>();
    if (!(a > T() && b > T() && is_int(a) && is_int(b) && a + b <= T(200))) return nn<T>();
    pc c = c_ibeta(static_cast<unsigned>(a), static_cast<unsigned>(b));
    if (comp) {
        for (auto& v : c) v = -v;
        c[0] += 1;
    }
    const std::complex<T> v = poly_ray_sym(c, x, y);
    return hit(v.real(), v.imag());
}
// zero of an odd polynomial at x = +-0: the sign of x times the sign of the linear coefficient
template <class T>
[[nodiscard]] inline quaternion<T> odd_zero(const quaternion<T>& q, long double c1) {
    const bool neg = std::signbit(q.w()) != (c1 < 0);
    return {neg ? -T() : T(), q.x(), q.y(), q.z()};
}
}  // namespace detail::sf::sv

// Polynomials: non-finite arguments take the limit along the ray (README), and an odd
// polynomial gives a signed zero at +-0.
#define QUAT_SF_POLY_SV(COEFFS, ODD)                                                          \
    if (quat::detail::needs_sv(q)) {                                                          \
        if (!quat::detail::is_finite(q))                                                      \
            return detail::sf::sv::poly_sv(q, true, [&](T x_, T y_) {                         \
                return detail::sf::sv::poly_ray_sym<T>(COEFFS, x_, y_);                       \
            });                                                                               \
        if (q.w() == T() && (ODD)) return detail::sf::sv::odd_zero(q, (COEFFS)[1]);           \
    }
/// Legendre polynomial P_l(q).
template <class T>
[[nodiscard]] inline quaternion<T> legendre(unsigned l, const quaternion<T>& q) {
    QUAT_SF_W;
    QUAT_SF_POLY_SV(detail::sf::sv::c_legendre(l), l % 2 == 1)
    const auto x = detail::sf::to_w(q);
    quaternion<W> p0(W(1)), p1 = x;
    if (l == 0) return detail::sf::from_w<T>(p0);
    for (unsigned k = 1; k < l; ++k) {
        quaternion<W> p2 = (W(2 * k + 1) * x * p1 - W(k) * p0) / W(k + 1);
        p0 = p1;
        p1 = p2;
    }
    return detail::sf::from_w<T>(p1);
}
/// Associated Legendre function P_l^m(q) = (1 - q^2)^{m/2} d^m/dq^m P_l(q), *without* the
/// Condon-Shortley phase (the <cmath> std::assoc_legendre convention).  For odd m the
/// principal square root of 1 - q^2 is used.
template <class T>
[[nodiscard]] inline quaternion<T> assoc_legendre(unsigned l, unsigned m, const quaternion<T>& q) {
    QUAT_SF_W;
    if (m > l) return quaternion<T>();
    if (quat::detail::needs_sv(q)) {
        using detail::sf::sv::pc;
        const pc d = detail::sf::sv::poly_deriv(detail::sf::sv::c_legendre(l), m);
        if (!quat::detail::is_finite(q)) {
            if (m % 2 == 0) {
                const pc c = detail::sf::sv::poly_mul(detail::sf::sv::c_one_minus_x2(m / 2), d);
                return detail::sf::sv::poly_sv(q, true, [&](T x_, T y_) { return detail::sf::sv::poly_ray_sym<T>(c, x_, y_); });
            }
            // odd m: P = s R(z), s = sqrt(1 - z^2) ~ -i z sqrt(1 - z^-2) for Im z > 0; the
            // limit is that of the polynomial part of z sqrt(1 - z^-2) R(z)
            const pc r = detail::sf::sv::poly_mul(detail::sf::sv::c_one_minus_x2((m - 1) / 2), d);
            pc t(r.size() + 1, 0);
            for (std::size_t p = 0; p < t.size(); ++p) {
                long double b = 1;  // coefficients of sqrt(1 - w)
                for (std::size_t j = 0; p + 2 * j <= r.size(); ++j) {
                    if (p + 2 * j >= 1) t[p] += b * r[p + 2 * j - 1];
                    b *= (static_cast<long double>(j) - 0.5L) / static_cast<long double>(j + 1);
                }
            }
            return detail::sf::sv::poly_sv(q, false, [&](T x_, T y_) -> std::complex<T> {
                // on the real axis (|x| > 1) s = +i sqrt(x^2 - 1), as for finite x
                if (y_ == T()) return {T(), (x_ < T() ? T(-1) : T(1)) * detail::sf::sv::poly_ray<T>(t, x_, T()).real()};
                const std::complex<T> v = detail::sf::sv::poly_ray<T>(t, x_, std::abs(y_));
                const std::complex<T> o(v.imag(), -v.real());
                return std::signbit(y_) ? std::conj(o) : o;
            });
        }
        if (q.w() == T() && (l - m) % 2 == 1) return detail::sf::sv::odd_zero(q, d[1]);
    }
    const auto x = detail::sf::to_w(q);
    const quaternion<W> s = sqrt(W(1) - x * x);
    // P_m^m = (2m-1)!! s^m
    quaternion<W> pmm(W(1));
    for (unsigned i = 1; i <= m; ++i) pmm = pmm * s * W(2 * i - 1);
    if (l == m) return detail::sf::from_w<T>(pmm);
    quaternion<W> pm1 = W(2 * m + 1) * x * pmm;
    for (unsigned k = m + 2; k <= l; ++k) {
        quaternion<W> p2 = (W(2 * k - 1) * x * pm1 - W(k + m - 1) * pmm) / W(k - m);
        pmm = pm1;
        pm1 = p2;
    }
    return detail::sf::from_w<T>(pm1);
}
/// Spherical harmonic angular part Y_l^m(theta, 0) for a quaternion polar angle theta
/// (the <cmath> std::sph_legendre convention, including the Condon-Shortley phase).
template <class T>
[[nodiscard]] inline quaternion<T> sph_legendre(unsigned l, unsigned m, const quaternion<T>& theta) {
    QUAT_SF_W;
    if (m > l) return quaternion<T>();
    if (quat::detail::needs_sv(theta)) {
        if (!quat::detail::is_finite(theta)) {
            const T c00 = T(0.5) / detail::lm::sqrt(quat::detail::pi_v<T>);
            return detail::sf::sv::poly_sv(theta, true, [&](T x_, T y_) -> std::complex<T> {
                if (l == 0) return {c00, detail::sf::sv::zs(y_)};
                if (!std::isinf(y_)) return {detail::sf::nan<T>(), y_ == T() ? y_ : detail::sf::nan<T>()};
                if (std::isinf(x_)) return {detail::sf::nan<T>(), detail::sf::nan<T>()};
                // Y ~ N (-1)^m (i e^{-i theta} / 2)^m (e^{-i theta} / 2)^{l-m} P_l^{(m)} leading term
                const std::complex<T> c = detail::sf::sv::ph(T(m) / 2, T(-1), T(l) * x_);
                const T sg = m % 2 ? T(-1) : T(1);
                return {detail::sf::sv::ci(sg * c.real()), detail::sf::sv::ci(sg * c.imag())};
            });
        }
        if (theta.w() == T() && m > 0) {
            const bool neg = (m % 2 == 1) && !std::signbit(theta.w());
            return {neg ? -T() : T(), theta.x(), theta.y(), theta.z()};
        }
    }
    const auto t = detail::sf::to_w(theta);
    const quaternion<W> x = cos(t), s = sin(t);
    // normalized recurrence: Ybar_l^m = N_l^m P_l^m with N = sqrt((2l+1)/(4pi) (l-m)!/(l+m)!)
    const W pi = detail::pi_v<W>;
    W lnorm = W(0.5) * (std::log(W(2 * m + 1) / (W(4) * pi)) - detail::lm::lgamma(W(2 * m + 1)));
    // (2m-1)!! in log form
    for (unsigned i = 1; i <= m; ++i) lnorm += std::log(W(2 * i - 1));
    quaternion<W> ymm = detail::lm::exp(lnorm) * pow(s, static_cast<int>(m));
    if (m % 2 == 1) ymm = -ymm;
    if (l == m) return detail::sf::from_w<T>(ymm);
    quaternion<W> y1 = detail::lm::sqrt(W(2 * m + 3)) * x * ymm;
    for (unsigned k = m + 2; k <= l; ++k) {
        const W kk = W(k), mm = W(m);
        const W a = detail::lm::sqrt((W(4) * kk * kk - W(1)) / (kk * kk - mm * mm));
        const W b = detail::lm::sqrt(((kk - W(1)) * (kk - W(1)) - mm * mm) / (W(4) * (kk - W(1)) * (kk - W(1)) - W(1)));
        quaternion<W> y2 = a * (x * y1 - b * ymm);
        ymm = y1;
        y1 = y2;
    }
    return detail::sf::from_w<T>(y1);
}
/// Physicists' Hermite polynomial H_n(q).
template <class T>
[[nodiscard]] inline quaternion<T> hermite(unsigned n, const quaternion<T>& q) {
    QUAT_SF_W;
    QUAT_SF_POLY_SV(detail::sf::sv::c_hermite(n), n % 2 == 1)
    const auto x = detail::sf::to_w(q);
    quaternion<W> h0(W(1)), h1 = W(2) * x;
    if (n == 0) return detail::sf::from_w<T>(h0);
    for (unsigned k = 1; k < n; ++k) {
        quaternion<W> h2 = W(2) * x * h1 - W(2 * k) * h0;
        h0 = h1;
        h1 = h2;
    }
    return detail::sf::from_w<T>(h1);
}
/// Associated (generalized) Laguerre polynomial L_n^m(q); laguerre(n, q) = L_n^0(q).
template <class T>
[[nodiscard]] inline quaternion<T> assoc_laguerre(unsigned n, unsigned m, const quaternion<T>& q) {
    QUAT_SF_W;
    QUAT_SF_POLY_SV(detail::sf::sv::c_laguerre(n, m), false)
    const auto x = detail::sf::to_w(q);
    const W a = W(m);
    quaternion<W> l0(W(1)), l1 = W(1) + a - x;
    if (n == 0) return detail::sf::from_w<T>(l0);
    for (unsigned k = 1; k < n; ++k) {
        quaternion<W> l2 = ((W(2 * k + 1) + a - x) * l1 - (W(k) + a) * l0) / W(k + 1);
        l0 = l1;
        l1 = l2;
    }
    return detail::sf::from_w<T>(l1);
}
template <class T>
[[nodiscard]] inline quaternion<T> laguerre(unsigned n, const quaternion<T>& q) {
    return assoc_laguerre(n, 0u, q);
}
/// Chebyshev polynomials of the first and second kind T_n(q), U_n(q).
template <class T>
[[nodiscard]] inline quaternion<T> chebyshev_t(unsigned n, const quaternion<T>& q) {
    QUAT_SF_W;
    QUAT_SF_POLY_SV(detail::sf::sv::c_chebyshev(n, false), n % 2 == 1)
    const auto x = detail::sf::to_w(q);
    quaternion<W> t0(W(1)), t1 = x;
    if (n == 0) return detail::sf::from_w<T>(t0);
    for (unsigned k = 1; k < n; ++k) {
        quaternion<W> t2 = W(2) * x * t1 - t0;
        t0 = t1;
        t1 = t2;
    }
    return detail::sf::from_w<T>(t1);
}
template <class T>
[[nodiscard]] inline quaternion<T> chebyshev_u(unsigned n, const quaternion<T>& q) {
    QUAT_SF_W;
    QUAT_SF_POLY_SV(detail::sf::sv::c_chebyshev(n, true), n % 2 == 1)
    const auto x = detail::sf::to_w(q);
    quaternion<W> u0(W(1)), u1 = W(2) * x;
    if (n == 0) return detail::sf::from_w<T>(u0);
    for (unsigned k = 1; k < n; ++k) {
        quaternion<W> u2 = W(2) * x * u1 - u0;
        u0 = u1;
        u1 = u2;
    }
    return detail::sf::from_w<T>(u1);
}

// ---- elliptic integrals (std convention: modulus k, characteristic nu in 1 - nu sin^2) ----
/// K(k) for a quaternion modulus k.
template <class T>
[[nodiscard]] inline quaternion<T> comp_ellint_1(const quaternion<T>& k) {
    QUAT_SF_W;
    return detail::sf::lift_s(k, [](const QUAT_SF_C& z) { return detail::sf::comp_ellint_1k(z); },
                            [](W x) { return std::abs(x) <= W(1); },
        [](T x, T y) { return detail::sf::sv::r_to_zero(x, y); });
}
/// E(k) for a quaternion modulus k.
template <class T>
[[nodiscard]] inline quaternion<T> comp_ellint_2(const quaternion<T>& k) {
    QUAT_SF_W;
    return detail::sf::lift_s(k, [](const QUAT_SF_C& z) { return detail::sf::comp_ellint_2k(z); },
                            [](W x) { return std::abs(x) <= W(1); },
        [](T x, T y) { return detail::sf::sv::r_comp_ellint_2(x, y); });
}
/// Pi(nu, k) for real modulus |k| <= 1 and a quaternion characteristic nu.
template <class T>
[[nodiscard]] inline quaternion<T> comp_ellint_3(std::type_identity_t<T> k, const quaternion<T>& nu) {
    QUAT_SF_W;
    const W kk = static_cast<W>(k);
    return detail::sf::lift_s(nu, [kk](const QUAT_SF_C& z) { return detail::sf::comp_ellint_3n(kk, z); },
                            [](W x) { return x < W(1); },
        [](T x, T y) { return detail::sf::sv::r_to_zero(x, y); });
}
/// F(k, phi) for real |k| <= 1 and a quaternion amplitude phi.
template <class T>
[[nodiscard]] inline quaternion<T> ellint_1(std::type_identity_t<T> k, const quaternion<T>& phi) {
    QUAT_SF_W;
    const W kk = static_cast<W>(k);
    return detail::sf::lift_s(phi, [kk](const QUAT_SF_C& z) { return detail::sf::ellint_1p(kk, z); },
                            [kk](W) { return std::abs(kk) <= W(1); },
        [k](T x, T y) { return detail::sf::sv::r_ellint_1(static_cast<T>(k), x, y); });
}
/// E(k, phi) for real |k| <= 1 and a quaternion amplitude phi.
template <class T>
[[nodiscard]] inline quaternion<T> ellint_2(std::type_identity_t<T> k, const quaternion<T>& phi) {
    QUAT_SF_W;
    const W kk = static_cast<W>(k);
    return detail::sf::lift_s(phi, [kk](const QUAT_SF_C& z) { return detail::sf::ellint_2p(kk, z); },
                            [kk](W) { return std::abs(kk) <= W(1); },
        [k](T x, T y) { return detail::sf::sv::r_ellint_2(static_cast<T>(k), x, y); });
}
/// Pi(nu, k, phi) for real k, real phi and a quaternion characteristic nu.
template <class T>
[[nodiscard]] inline quaternion<T> ellint_3(std::type_identity_t<T> k, const quaternion<T>& nu,
                                            std::type_identity_t<T> phi) {
    QUAT_SF_W;
    const W kk = static_cast<W>(k), ph = static_cast<W>(phi);
    return detail::sf::lift_s(nu, [kk, ph](const QUAT_SF_C& z) { return detail::sf::ellint_3n(kk, z, ph); },
                            [ph](W x) { const W s = detail::lm::sin(ph); return x * s * s < W(1); },
        [](T x, T y) { return detail::sf::sv::r_to_zero(x, y); });
}

// ---- exponential, sine/cosine and Fresnel integrals ----------------------------------------
/// Ei(q) (exponential integral; on the real axis the real principal value as std::expint).
template <class T>
[[nodiscard]] inline quaternion<T> expint(const quaternion<T>& q) {
    QUAT_SF_W;
    return detail::sf::lift_s(q, [](const QUAT_SF_C& z) { return detail::sf::expint_ei(z); }, true,
        [](T x, T y) { return detail::sf::sv::r_ei(x, y); });
}
/// E_n(q) = int_1^inf e^{-qt} t^{-n} dt (E_1 is the exponential integral E1).
template <class T>
[[nodiscard]] inline quaternion<T> expint(unsigned n, const quaternion<T>& q) {
    QUAT_SF_W;
    return detail::sf::lift_s(
        q,
        [n](const QUAT_SF_C& z) -> QUAT_SF_C {
            using namespace detail::sf;
            if (n == 0) return detail::lm::exp(-z) / z;
            if (n == 1) return expint_e1(z);
            if (z == QUAT_SF_C()) return {W(1) / W(n - 1), W()};
            return detail::lm::exp(W(n - 1) * std::log(z)) * gamma_upper(W(1) - W(n), z);
        },
        [n](W x) { return x > W() || (n > 1 && x == W()); },
        [n](T x, T y) { return detail::sf::sv::r_en(n, x, y); });
}
/// Sine integral Si(q) (entire) and cosine integral Ci(q).
template <class T>
[[nodiscard]] inline quaternion<T> sinint(const quaternion<T>& q) {
    QUAT_SF_W;
    return detail::sf::lift_s(q, [](const QUAT_SF_C& z) { return detail::sf::sici(z, false); }, true,
        [](T x, T y) { return detail::sf::sv::r_si(x, y); });
}
template <class T>
[[nodiscard]] inline quaternion<T> cosint(const quaternion<T>& q) {
    QUAT_SF_W;
    return detail::sf::lift_s(q, [](const QUAT_SF_C& z) { return detail::sf::sici(z, true); },
                            [](W x) { return x > W(); },
        [](T x, T y) { return detail::sf::sv::r_ci(x, y); });
}
/// Fresnel integrals S(q) = int_0^q sin(pi t^2/2) dt and C(q) (entire).
template <class T>
[[nodiscard]] inline quaternion<T> fresnel_s(const quaternion<T>& q) {
    QUAT_SF_W;
    return detail::sf::lift_s(q, [](const QUAT_SF_C& z) { return detail::sf::fresnel(z, true); }, true,
        [](T x, T y) { return detail::sf::sv::r_fresnel(true, x, y); });
}
template <class T>
[[nodiscard]] inline quaternion<T> fresnel_c(const quaternion<T>& q) {
    QUAT_SF_W;
    return detail::sf::lift_s(q, [](const QUAT_SF_C& z) { return detail::sf::fresnel(z, false); }, true,
        [](T x, T y) { return detail::sf::sv::r_fresnel(false, x, y); });
}

// ---- polygamma, Hurwitz zeta, polylogarithm, Lambert W ------------------------------------
/// psi^(m)(q), the m-th derivative of the digamma function.
template <class T>
[[nodiscard]] inline quaternion<T> polygamma(unsigned m, const quaternion<T>& q) {
    QUAT_SF_W;
    return detail::sf::lift_s(q, [m](const QUAT_SF_C& z) { return detail::sf::polygamma_k(m, z); }, true,
        [m](T x, T y) { return detail::sf::sv::r_polygamma(m, x, y); });
}
template <class T>
[[nodiscard]] inline quaternion<T> trigamma(const quaternion<T>& q) {
    return polygamma(1u, q);
}
/// Hurwitz zeta zeta(s, a) for a quaternion s and real a (a not a non-positive integer).
template <class T>
[[nodiscard]] inline quaternion<T> hurwitz_zeta(const quaternion<T>& s, std::type_identity_t<T> a) {
    QUAT_SF_W;
    const W wa = static_cast<W>(a);
    return detail::sf::lift_s(s, [wa](const QUAT_SF_C& z) { return detail::sf::hurwitz_zeta_k(z, QUAT_SF_C(wa)); },
                            [wa](W) { return wa > W(); },
        [a](T x, T y) { return detail::sf::sv::r_hurwitz(static_cast<T>(a), x, y); });
}
/// Polylogarithm Li_s(q) for real s (cut [1, inf) on the real axis).
template <class T, class N>
    requires std::is_arithmetic_v<N>
[[nodiscard]] inline quaternion<T> polylog(N s, const quaternion<T>& q) {
    QUAT_SF_W;
    const W ws = static_cast<W>(s);
    return detail::sf::lift_s(q, [ws](const QUAT_SF_C& z) { return detail::sf::polylog_k(ws, z); },
                            [](W x) { return x <= W(1); },
        [s](T x, T y) { return detail::sf::sv::r_polylog(static_cast<T>(s), x, y); });
}
/// Lambert W, branch k (k = 0 is the principal branch).
template <class T>
[[nodiscard]] inline quaternion<T> lambert_w(int k, const quaternion<T>& q) {
    QUAT_SF_W;
    return detail::sf::lift_s(q, [k](const QUAT_SF_C& z) { return detail::sf::lambert_w_k(k, z); },
                            [k](W x) {
                                const W e1 = W(1) / detail::lm::exp(W(1));
                                if (k == 0) return x >= -e1;
                                if (k == -1) return x >= -e1 && x < W();
                                return false;
                            },
        [k](T x, T y) { return detail::sf::sv::r_lambert_w(k, x, y); },
        [k](W x) { return std::conj(detail::sf::lambert_w_k(-k, QUAT_SF_C(x, W()))); });
}
template <class T>
[[nodiscard]] inline quaternion<T> lambert_w(const quaternion<T>& q) {
    return lambert_w(0, q);
}

// ---- incomplete gamma and beta (real parameters) ------------------------------------------
/// Upper incomplete gamma Gamma(a, q) and lower gamma(a, q), real a.
template <class T>
[[nodiscard]] inline quaternion<T> tgamma_upper(std::type_identity_t<T> a, const quaternion<T>& q) {
    QUAT_SF_W;
    const W wa = static_cast<W>(a);
    return detail::sf::lift_s(q, [wa](const QUAT_SF_C& z) { return detail::sf::gamma_upper(wa, z); },
                            [wa](W x) { return x > W() || (detail::sf::is_integer(wa) && wa > W()); },
        [a](T x, T y) { return detail::sf::sv::r_gamma_inc(static_cast<T>(a), true, false, x, y); });
}
template <class T>
[[nodiscard]] inline quaternion<T> tgamma_lower(std::type_identity_t<T> a, const quaternion<T>& q) {
    QUAT_SF_W;
    const W wa = static_cast<W>(a);
    return detail::sf::lift_s(q, [wa](const QUAT_SF_C& z) { return detail::sf::gamma_lower(wa, z); },
                            [wa](W x) { return x >= W() || (detail::sf::is_integer(wa) && wa > W()); },
        [a](T x, T y) { return detail::sf::sv::r_gamma_inc(static_cast<T>(a), false, false, x, y); });
}
/// Regularized P(a, q) = gamma(a, q)/Gamma(a) and Q(a, q) = Gamma(a, q)/Gamma(a), a > 0.
template <class T>
[[nodiscard]] inline quaternion<T> gamma_p(std::type_identity_t<T> a, const quaternion<T>& q) {
    QUAT_SF_W;
    const W wa = static_cast<W>(a);
    return detail::sf::lift_s(
        q,
        [wa](const QUAT_SF_C& z) {
            return wa > W() ? detail::sf::gamma_pq(wa, z).p : QUAT_SF_C(detail::sf::nan<W>(), W());
        },
        [wa](W x) { return x >= W() || detail::sf::is_integer(wa); },
        [a](T x, T y) { return detail::sf::sv::r_gamma_inc(static_cast<T>(a), false, true, x, y); });
}
template <class T>
[[nodiscard]] inline quaternion<T> gamma_q(std::type_identity_t<T> a, const quaternion<T>& q) {
    QUAT_SF_W;
    const W wa = static_cast<W>(a);
    return detail::sf::lift_s(
        q,
        [wa](const QUAT_SF_C& z) {
            return wa > W() ? detail::sf::gamma_pq(wa, z).q : QUAT_SF_C(detail::sf::nan<W>(), W());
        },
        [wa](W x) { return x >= W() || detail::sf::is_integer(wa); },
        [a](T x, T y) { return detail::sf::sv::r_gamma_inc(static_cast<T>(a), true, true, x, y); });
}
/// Regularized incomplete beta I_q(a, b) and its complement 1 - I_q(a, b), real a, b > 0.
template <class T>
[[nodiscard]] inline quaternion<T> ibeta(std::type_identity_t<T> a, std::type_identity_t<T> b,
                                         const quaternion<T>& q) {
    QUAT_SF_W;
    const W wa = static_cast<W>(a), wb = static_cast<W>(b);
    return detail::sf::lift_s(q, [wa, wb](const QUAT_SF_C& z) { return detail::sf::ibeta_k(wa, wb, z); },
                            [](W x) { return x >= W() && x <= W(1); },
        [a, b](T x, T y) { return detail::sf::sv::r_ibeta(static_cast<T>(a), static_cast<T>(b), false, x, y); });
}
template <class T>
[[nodiscard]] inline quaternion<T> ibetac(std::type_identity_t<T> a, std::type_identity_t<T> b,
                                          const quaternion<T>& q) {
    QUAT_SF_W;
    const W wa = static_cast<W>(a), wb = static_cast<W>(b);
    return detail::sf::lift_s(
        q, [wa, wb](const QUAT_SF_C& z) { return W(1) - detail::sf::ibeta_k(wa, wb, z); },
        [](W x) { return x >= W() && x <= W(1); },
        [a, b](T x, T y) { return detail::sf::sv::r_ibeta(static_cast<T>(a), static_cast<T>(b), true, x, y); });
}

#undef QUAT_SF_W
#undef QUAT_SF_C

}  // namespace quat

#if defined(QUATERNION_MSVC_CONTRACT_OFF) && !defined(QUATERNION_DIAG_NORESTORE)
#pragma float_control(pop)
#endif

#endif  // QUATERNION_SPECIAL_FUNCTIONS_HPP
