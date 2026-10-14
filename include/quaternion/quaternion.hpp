// quaternion.hpp - header-only C++23 quaternion numbers modeled on <complex>.
// SPDX-License-Identifier: MIT
//
// A quaternion q = a + b*i + c*j + d*k, with i^2 = j^2 = k^2 = ijk = -1.
// Multiplication is associative but NOT commutative.
#ifndef QUATERNION_QUATERNION_HPP
#define QUATERNION_QUATERNION_HPP

#include <cmath>
#include <complex>
#include <concepts>
#include <cstddef>
#include <format>
#include <istream>
#include <limits>
#include <numbers>
#include <ostream>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

#define QUATERNION_VERSION_MAJOR 1
#define QUATERNION_VERSION_MINOR 0
#define QUATERNION_VERSION_PATCH 0

namespace quat {

template <std::floating_point T>
class quaternion;

namespace detail {
template <class>
inline constexpr bool is_quaternion_v = false;
template <class T>
inline constexpr bool is_quaternion_v<quaternion<T>> = true;
}  // namespace detail

/// Quaternion number a + bi + cj + dk over the floating-point type T.
/// Like std::complex, only float, double and long double are supported.
template <std::floating_point T>
class quaternion {
public:
    using value_type = T;

    // ---- construction ---------------------------------------------------
    constexpr quaternion(const T& a = T(), const T& b = T(), const T& c = T(),
                         const T& d = T()) noexcept
        : a_(a), b_(b), c_(c), d_(d) {}

    /// From a pair of complex numbers: z1 + z2*j  ==  (z1.re, z1.im, z2.re, z2.im).
    constexpr quaternion(const std::complex<T>& z1,
                         const std::complex<T>& z2 = std::complex<T>()) noexcept
        : a_(z1.real()), b_(z1.imag()), c_(z2.real()), d_(z2.imag()) {}

    constexpr quaternion(const quaternion&) noexcept = default;
    constexpr quaternion& operator=(const quaternion&) noexcept = default;

    /// Converting constructor; explicit when the conversion narrows (as std::complex).
    template <std::floating_point U>
    constexpr explicit(!requires(U u) { T{u}; })
        quaternion(const quaternion<U>& q) noexcept
        : a_(static_cast<T>(q.R_component_1())), b_(static_cast<T>(q.R_component_2())),
          c_(static_cast<T>(q.R_component_3())), d_(static_cast<T>(q.R_component_4())) {}

    // ---- accessors ------------------------------------------------------
    [[nodiscard]] constexpr T real() const noexcept { return a_; }
    constexpr void real(T v) noexcept { a_ = v; }

    /// Vector (unreal) part as a pure quaternion (0, b, c, d).
    [[nodiscard]] constexpr quaternion unreal() const noexcept { return {T(), b_, c_, d_}; }

    [[nodiscard]] constexpr T R_component_1() const noexcept { return a_; }
    [[nodiscard]] constexpr T R_component_2() const noexcept { return b_; }
    [[nodiscard]] constexpr T R_component_3() const noexcept { return c_; }
    [[nodiscard]] constexpr T R_component_4() const noexcept { return d_; }
    [[nodiscard]] constexpr std::complex<T> C_component_1() const noexcept { return {a_, b_}; }
    [[nodiscard]] constexpr std::complex<T> C_component_2() const noexcept { return {c_, d_}; }

    [[nodiscard]] constexpr T w() const noexcept { return a_; }
    [[nodiscard]] constexpr T x() const noexcept { return b_; }
    [[nodiscard]] constexpr T y() const noexcept { return c_; }
    [[nodiscard]] constexpr T z() const noexcept { return d_; }
    constexpr void w(T v) noexcept { a_ = v; }
    constexpr void x(T v) noexcept { b_ = v; }
    constexpr void y(T v) noexcept { c_ = v; }
    constexpr void z(T v) noexcept { d_ = v; }

    /// Tuple-like element access (also used by structured bindings).
    template <std::size_t I>
        requires(I < 4)
    [[nodiscard]] constexpr T& get() & noexcept { return comp<I>(*this); }
    template <std::size_t I>
        requires(I < 4)
    [[nodiscard]] constexpr const T& get() const& noexcept { return comp<I>(*this); }
    template <std::size_t I>
        requires(I < 4)
    [[nodiscard]] constexpr T&& get() && noexcept { return std::move(comp<I>(*this)); }
    template <std::size_t I>
        requires(I < 4)
    [[nodiscard]] constexpr const T&& get() const&& noexcept { return std::move(comp<I>(*this)); }

    // ---- assignment -----------------------------------------------------
    constexpr quaternion& operator=(const T& s) noexcept {
        a_ = s; b_ = c_ = d_ = T();
        return *this;
    }
    constexpr quaternion& operator=(const std::complex<T>& z) noexcept {
        a_ = z.real(); b_ = z.imag(); c_ = d_ = T();
        return *this;
    }
    template <std::floating_point U>
    constexpr quaternion& operator=(const quaternion<U>& q) noexcept {
        a_ = static_cast<T>(q.R_component_1()); b_ = static_cast<T>(q.R_component_2());
        c_ = static_cast<T>(q.R_component_3()); d_ = static_cast<T>(q.R_component_4());
        return *this;
    }

    // ---- compound assignment: scalar ---------------------------------------
    constexpr quaternion& operator+=(const T& s) noexcept { a_ += s; return *this; }
    constexpr quaternion& operator-=(const T& s) noexcept { a_ -= s; return *this; }
    constexpr quaternion& operator*=(const T& s) noexcept {
        a_ *= s; b_ *= s; c_ *= s; d_ *= s;
        return *this;
    }
    constexpr quaternion& operator/=(const T& s) noexcept {
        a_ /= s; b_ /= s; c_ /= s; d_ /= s;
        return *this;
    }

    // ---- compound assignment: complex (treated as a + bi) --------------------
    constexpr quaternion& operator+=(const std::complex<T>& z) noexcept { return *this += quaternion(z); }
    constexpr quaternion& operator-=(const std::complex<T>& z) noexcept { return *this -= quaternion(z); }
    constexpr quaternion& operator*=(const std::complex<T>& z) noexcept { return *this *= quaternion(z); }
    constexpr quaternion& operator/=(const std::complex<T>& z) noexcept { return *this /= quaternion(z); }

    // ---- compound assignment: quaternion ----------------------------------
    template <std::floating_point U>
    constexpr quaternion& operator+=(const quaternion<U>& q) noexcept {
        a_ += static_cast<T>(q.R_component_1()); b_ += static_cast<T>(q.R_component_2());
        c_ += static_cast<T>(q.R_component_3()); d_ += static_cast<T>(q.R_component_4());
        return *this;
    }
    template <std::floating_point U>
    constexpr quaternion& operator-=(const quaternion<U>& q) noexcept {
        a_ -= static_cast<T>(q.R_component_1()); b_ -= static_cast<T>(q.R_component_2());
        c_ -= static_cast<T>(q.R_component_3()); d_ -= static_cast<T>(q.R_component_4());
        return *this;
    }
    /// Right-multiplication: *this = *this * q (Hamilton product).
    template <std::floating_point U>
    constexpr quaternion& operator*=(const quaternion<U>& q) noexcept {
        const T a2 = static_cast<T>(q.R_component_1()), b2 = static_cast<T>(q.R_component_2());
        const T c2 = static_cast<T>(q.R_component_3()), d2 = static_cast<T>(q.R_component_4());
        const T a = a_ * a2 - b_ * b2 - c_ * c2 - d_ * d2;
        const T b = a_ * b2 + b_ * a2 + c_ * d2 - d_ * c2;
        const T c = a_ * c2 - b_ * d2 + c_ * a2 + d_ * b2;
        const T d = a_ * d2 + b_ * c2 - c_ * b2 + d_ * a2;
        a_ = a; b_ = b; c_ = c; d_ = d;
        return *this;
    }
    /// Right-division: *this = *this * q^-1.
    template <std::floating_point U>
    constexpr quaternion& operator/=(const quaternion<U>& q) noexcept {
        const quaternion r(q);
        const T n = r.a_ * r.a_ + r.b_ * r.b_ + r.c_ * r.c_ + r.d_ * r.d_;
        *this *= quaternion(r.a_ / n, -r.b_ / n, -r.c_ / n, -r.d_ / n);
        return *this;
    }

    // ---- unary ------------------------------------------------------------
    [[nodiscard]] constexpr quaternion operator+() const noexcept { return *this; }
    [[nodiscard]] constexpr quaternion operator-() const noexcept { return {-a_, -b_, -c_, -d_}; }

    // ---- binary (hidden friends, so implicit conversions such as int -> T apply)
#define QUAT_BINOP(OP)                                                                      \
    [[nodiscard]] friend constexpr quaternion operator OP(const quaternion& l,              \
                                                          const quaternion& r) noexcept {   \
        quaternion t(l); t OP## = r; return t;                                             \
    }                                                                                       \
    [[nodiscard]] friend constexpr quaternion operator OP(const quaternion& l, const T& r) noexcept { \
        quaternion t(l); t OP## = r; return t;                                             \
    }                                                                                       \
    [[nodiscard]] friend constexpr quaternion operator OP(const T& l, const quaternion& r) noexcept { \
        quaternion t(l); t OP## = r; return t;                                             \
    }                                                                                       \
    [[nodiscard]] friend constexpr quaternion operator OP(const quaternion& l,              \
                                                          const std::complex<T>& r) noexcept { \
        quaternion t(l); t OP## = r; return t;                                             \
    }                                                                                       \
    [[nodiscard]] friend constexpr quaternion operator OP(const std::complex<T>& l,         \
                                                          const quaternion& r) noexcept {   \
        quaternion t(l); t OP## = r; return t;                                             \
    }
    QUAT_BINOP(+)
    QUAT_BINOP(-)
    QUAT_BINOP(*)
    QUAT_BINOP(/)
#undef QUAT_BINOP

    // ---- comparison (C++20 rewrites give !=, and reversed operands) ---------
    [[nodiscard]] friend constexpr bool operator==(const quaternion& l, const quaternion& r) noexcept {
        return l.a_ == r.a_ && l.b_ == r.b_ && l.c_ == r.c_ && l.d_ == r.d_;
    }
    [[nodiscard]] friend constexpr bool operator==(const quaternion& l, const T& r) noexcept {
        return l.a_ == r && l.b_ == T() && l.c_ == T() && l.d_ == T();
    }
    [[nodiscard]] friend constexpr bool operator==(const quaternion& l, const std::complex<T>& r) noexcept {
        return l.a_ == r.real() && l.b_ == r.imag() && l.c_ == T() && l.d_ == T();
    }

private:
    template <std::size_t I, class Self>
    static constexpr auto& comp(Self& s) noexcept {
        if constexpr (I == 0) return s.a_;
        else if constexpr (I == 1) return s.b_;
        else if constexpr (I == 2) return s.c_;
        else return s.d_;
    }

    T a_, b_, c_, d_;
};

// CTAD: quaternion(1.0, 2.0, 3.0, 4.0) -> quaternion<double>
template <std::floating_point T>
quaternion(T) -> quaternion<T>;
template <std::floating_point T>
quaternion(T, T) -> quaternion<T>;
template <std::floating_point T>
quaternion(T, T, T) -> quaternion<T>;
template <std::floating_point T>
quaternion(T, T, T, T) -> quaternion<T>;
template <std::floating_point T>
quaternion(std::complex<T>) -> quaternion<T>;
template <std::floating_point T>
quaternion(std::complex<T>, std::complex<T>) -> quaternion<T>;

// ---- tuple-like free get (found by ADL) ------------------------------------
template <std::size_t I, class T>
[[nodiscard]] constexpr T& get(quaternion<T>& q) noexcept { return q.template get<I>(); }
template <std::size_t I, class T>
[[nodiscard]] constexpr const T& get(const quaternion<T>& q) noexcept { return q.template get<I>(); }
template <std::size_t I, class T>
[[nodiscard]] constexpr T&& get(quaternion<T>&& q) noexcept { return std::move(q).template get<I>(); }
template <std::size_t I, class T>
[[nodiscard]] constexpr const T&& get(const quaternion<T>&& q) noexcept { return std::move(q).template get<I>(); }

// ---- value operations ------------------------------------------------------
template <class T>
[[nodiscard]] constexpr T real(const quaternion<T>& q) noexcept { return q.real(); }
template <class T>
[[nodiscard]] constexpr quaternion<T> unreal(const quaternion<T>& q) noexcept { return q.unreal(); }

/// Squared magnitude (field norm), a^2 + b^2 + c^2 + d^2, like std::norm.
template <class T>
[[nodiscard]] constexpr T norm(const quaternion<T>& q) noexcept {
    return q.w() * q.w() + q.x() * q.x() + q.y() * q.y() + q.z() * q.z();
}
template <class T>
[[nodiscard]] constexpr quaternion<T> conj(const quaternion<T>& q) noexcept {
    return {q.w(), -q.x(), -q.y(), -q.z()};
}
/// Multiplicative inverse conj(q)/norm(q).
template <class T>
[[nodiscard]] constexpr quaternion<T> inverse(const quaternion<T>& q) noexcept {
    return conj(q) / norm(q);
}
/// Euclidean (4D) inner product.
template <class T>
[[nodiscard]] constexpr T dot(const quaternion<T>& p, const quaternion<T>& q) noexcept {
    return p.w() * q.w() + p.x() * q.x() + p.y() * q.y() + p.z() * q.z();
}
/// Cross product of the vector parts; returns a pure quaternion.
template <class T>
[[nodiscard]] constexpr quaternion<T> cross(const quaternion<T>& p, const quaternion<T>& q) noexcept {
    return {T(), p.y() * q.z() - p.z() * q.y(), p.z() * q.x() - p.x() * q.z(),
            p.x() * q.y() - p.y() * q.x()};
}

namespace detail {
template <class T>
[[nodiscard]] inline T vabs(const quaternion<T>& q) noexcept {
    return std::hypot(q.x(), q.y(), q.z());
}
// Rebuild s + u*v where u is the unit vector of q (i if q has no vector part).
template <class T>
[[nodiscard]] inline quaternion<T> along(const quaternion<T>& q, T vlen, T s, T v) noexcept {
    if (vlen == T()) return {s, v, T(), T()};
    const T f = v / vlen;
    return {s, q.x() * f, q.y() * f, q.z() * f};
}
}  // namespace detail

/// Magnitude, computed without undue overflow/underflow.
template <class T>
[[nodiscard]] inline T abs(const quaternion<T>& q) noexcept {
    return std::hypot(std::hypot(q.w(), q.x()), std::hypot(q.y(), q.z()));
}
/// Angle in [0, pi] between q and the positive real axis: atan2(|v|, a).
template <class T>
[[nodiscard]] inline T arg(const quaternion<T>& q) noexcept {
    return std::atan2(detail::vabs(q), q.w());
}
/// Unit quaternion q/|q| (q itself if zero).
template <class T>
[[nodiscard]] inline quaternion<T> normalize(const quaternion<T>& q) noexcept {
    const T r = abs(q);
    return r == T() ? q : q / r;
}
/// Projection onto the Riemann-sphere analogue (as std::proj).
template <class T>
[[nodiscard]] inline quaternion<T> proj(const quaternion<T>& q) noexcept {
    if (std::isinf(q.w()) || std::isinf(q.x()) || std::isinf(q.y()) || std::isinf(q.z()))
        return {std::numeric_limits<T>::infinity(), std::copysign(T(), q.x()),
                std::copysign(T(), q.y()), std::copysign(T(), q.z())};
    return q;
}

/// rho * (cos(theta) + i sin(theta)) -- same as std::polar embedded on the i axis.
template <class T>
[[nodiscard]] inline quaternion<T> polar(const T& rho, const T& theta = T()) noexcept {
    return {rho * std::cos(theta), rho * std::sin(theta), T(), T()};
}
/// rho * (cos(theta) + u sin(theta)) with u = axis/|axis| (the vector part of axis).
template <class T>
[[nodiscard]] inline quaternion<T> polar(const T& rho, const T& theta,
                                         const quaternion<T>& axis) noexcept {
    return detail::along(axis, detail::vabs(axis), rho * std::cos(theta), rho * std::sin(theta));
}

// ---- transcendental functions ----------------------------------------------
template <class T>
[[nodiscard]] inline quaternion<T> exp(const quaternion<T>& q) noexcept {
    const T v = detail::vabs(q), e = std::exp(q.w());
    return detail::along(q, v, e * std::cos(v), e * std::sin(v));
}
/// Principal logarithm ln|q| + u*arg(q). For negative reals the result is ln|a| + pi*i.
template <class T>
[[nodiscard]] inline quaternion<T> log(const quaternion<T>& q) noexcept {
    const T v = detail::vabs(q);
    return detail::along(q, v, std::log(abs(q)), std::atan2(v, q.w()));
}
template <class T>
[[nodiscard]] inline quaternion<T> log10(const quaternion<T>& q) noexcept {
    return log(q) / std::numbers::ln10_v<T>;
}
/// Principal square root (non-negative real part).
template <class T>
[[nodiscard]] inline quaternion<T> sqrt(const quaternion<T>& q) noexcept {
    const T v = detail::vabs(q), r = abs(q);
    if (v == T()) {
        if (q.w() >= T()) return {std::sqrt(q.w()), q.x(), q.y(), q.z()};
        return {T(), std::sqrt(-q.w()), T(), T()};
    }
    const T s = std::sqrt((r + q.w()) / T(2));
    if (s == T()) return detail::along(q, v, T(), std::sqrt(r));
    return detail::along(q, v, s, v / (T(2) * s));
}

/// Integer power by repeated squaring (exact for small integer exponents).
template <class T>
[[nodiscard]] constexpr quaternion<T> pow(const quaternion<T>& q, int n) noexcept {
    if (n < 0) return inverse(pow(q, -(n + 1)) * q);
    quaternion<T> result(T(1)), base(q);
    for (unsigned k = static_cast<unsigned>(n); k; k >>= 1) {
        if (k & 1u) result *= base;
        base *= base;
    }
    return result;
}
template <class T>
[[nodiscard]] inline quaternion<T> pow(const quaternion<T>& q, const T& x) noexcept {
    if (q == quaternion<T>()) return x == T() ? quaternion<T>(T(1)) : quaternion<T>();
    return exp(log(q) * x);
}
/// q^p := exp(p * log(q)), the same operand order as std::pow for complex.
template <class T>
[[nodiscard]] inline quaternion<T> pow(const quaternion<T>& q, const quaternion<T>& p) noexcept {
    if (q == quaternion<T>()) return p == quaternion<T>() ? quaternion<T>(T(1)) : quaternion<T>();
    return exp(p * log(q));
}
template <class T>
[[nodiscard]] inline quaternion<T> pow(const T& x, const quaternion<T>& p) noexcept {
    if (x > T()) return exp(p * std::log(x));
    return pow(quaternion<T>(x), p);
}
template <class T>
[[nodiscard]] inline quaternion<T> pow(const quaternion<T>& q, const std::complex<T>& p) noexcept {
    return pow(q, quaternion<T>(p));
}
template <class T>
[[nodiscard]] inline quaternion<T> pow(const std::complex<T>& z, const quaternion<T>& p) noexcept {
    return pow(quaternion<T>(z), p);
}

template <class T>
[[nodiscard]] inline quaternion<T> sin(const quaternion<T>& q) noexcept {
    const T v = detail::vabs(q);
    return detail::along(q, v, std::sin(q.w()) * std::cosh(v), std::cos(q.w()) * std::sinh(v));
}
template <class T>
[[nodiscard]] inline quaternion<T> cos(const quaternion<T>& q) noexcept {
    const T v = detail::vabs(q);
    return detail::along(q, v, std::cos(q.w()) * std::cosh(v), -std::sin(q.w()) * std::sinh(v));
}
template <class T>
[[nodiscard]] inline quaternion<T> tan(const quaternion<T>& q) noexcept {
    return sin(q) * inverse(cos(q));  // sin(q) and cos(q) commute (same axis)
}
template <class T>
[[nodiscard]] inline quaternion<T> sinh(const quaternion<T>& q) noexcept {
    const T v = detail::vabs(q);
    return detail::along(q, v, std::sinh(q.w()) * std::cos(v), std::cosh(q.w()) * std::sin(v));
}
template <class T>
[[nodiscard]] inline quaternion<T> cosh(const quaternion<T>& q) noexcept {
    const T v = detail::vabs(q);
    return detail::along(q, v, std::cosh(q.w()) * std::cos(v), std::sinh(q.w()) * std::sin(v));
}
template <class T>
[[nodiscard]] inline quaternion<T> tanh(const quaternion<T>& q) noexcept {
    return sinh(q) * inverse(cosh(q));
}

// ---- stream I/O: "(a,b,c,d)" like std::complex -------------------------------
template <class T, class CharT, class Traits>
std::basic_ostream<CharT, Traits>& operator<<(std::basic_ostream<CharT, Traits>& os,
                                              const quaternion<T>& q) {
    std::basic_ostringstream<CharT, Traits> s;
    s.flags(os.flags());
    s.imbue(os.getloc());
    s.precision(os.precision());
    s << CharT('(') << q.w() << CharT(',') << q.x() << CharT(',') << q.y() << CharT(',')
      << q.z() << CharT(')');
    return os << s.str();
}

/// Accepts a, (a), (a,b), (a,b,c) and (a,b,c,d); sets failbit otherwise.
template <class T, class CharT, class Traits>
std::basic_istream<CharT, Traits>& operator>>(std::basic_istream<CharT, Traits>& is,
                                              quaternion<T>& q) {
    T v[4] = {T(), T(), T(), T()};
    CharT ch{};
    if (!(is >> ch)) return is;
    if (ch != CharT('(')) {
        is.putback(ch);
        if (is >> v[0]) q = quaternion<T>(v[0]);
        return is;
    }
    for (int n = 0; n < 4; ++n) {
        if (!(is >> v[n] >> ch)) return is;
        if (ch == CharT(')')) {
            q = quaternion<T>(v[0], v[1], v[2], v[3]);
            return is;
        }
        if (ch != CharT(',') || n == 3) break;
    }
    is.setstate(std::ios_base::failbit);
    return is;
}

// ---- user-defined literals ----------------------------------------------------
inline namespace literals {
inline namespace quaternion_literals {
[[nodiscard]] constexpr quaternion<double> operator""_i(long double v) noexcept {
    return {0.0, static_cast<double>(v), 0.0, 0.0};
}
[[nodiscard]] constexpr quaternion<double> operator""_i(unsigned long long v) noexcept {
    return {0.0, static_cast<double>(v), 0.0, 0.0};
}
[[nodiscard]] constexpr quaternion<double> operator""_j(long double v) noexcept {
    return {0.0, 0.0, static_cast<double>(v), 0.0};
}
[[nodiscard]] constexpr quaternion<double> operator""_j(unsigned long long v) noexcept {
    return {0.0, 0.0, static_cast<double>(v), 0.0};
}
[[nodiscard]] constexpr quaternion<double> operator""_k(long double v) noexcept {
    return {0.0, 0.0, 0.0, static_cast<double>(v)};
}
[[nodiscard]] constexpr quaternion<double> operator""_k(unsigned long long v) noexcept {
    return {0.0, 0.0, 0.0, static_cast<double>(v)};
}
[[nodiscard]] constexpr quaternion<float> operator""_if(long double v) noexcept {
    return {0.0f, static_cast<float>(v), 0.0f, 0.0f};
}
[[nodiscard]] constexpr quaternion<float> operator""_if(unsigned long long v) noexcept {
    return {0.0f, static_cast<float>(v), 0.0f, 0.0f};
}
[[nodiscard]] constexpr quaternion<float> operator""_jf(long double v) noexcept {
    return {0.0f, 0.0f, static_cast<float>(v), 0.0f};
}
[[nodiscard]] constexpr quaternion<float> operator""_jf(unsigned long long v) noexcept {
    return {0.0f, 0.0f, static_cast<float>(v), 0.0f};
}
[[nodiscard]] constexpr quaternion<float> operator""_kf(long double v) noexcept {
    return {0.0f, 0.0f, 0.0f, static_cast<float>(v)};
}
[[nodiscard]] constexpr quaternion<float> operator""_kf(unsigned long long v) noexcept {
    return {0.0f, 0.0f, 0.0f, static_cast<float>(v)};
}
[[nodiscard]] constexpr quaternion<long double> operator""_il(long double v) noexcept {
    return {0.0L, static_cast<long double>(v), 0.0L, 0.0L};
}
[[nodiscard]] constexpr quaternion<long double> operator""_il(unsigned long long v) noexcept {
    return {0.0L, static_cast<long double>(v), 0.0L, 0.0L};
}
[[nodiscard]] constexpr quaternion<long double> operator""_jl(long double v) noexcept {
    return {0.0L, 0.0L, static_cast<long double>(v), 0.0L};
}
[[nodiscard]] constexpr quaternion<long double> operator""_jl(unsigned long long v) noexcept {
    return {0.0L, 0.0L, static_cast<long double>(v), 0.0L};
}
[[nodiscard]] constexpr quaternion<long double> operator""_kl(long double v) noexcept {
    return {0.0L, 0.0L, 0.0L, static_cast<long double>(v)};
}
[[nodiscard]] constexpr quaternion<long double> operator""_kl(unsigned long long v) noexcept {
    return {0.0L, 0.0L, 0.0L, static_cast<long double>(v)};
}
}  // namespace quaternion_literals
}  // namespace literals

}  // namespace quat

// ---- tuple protocol -------------------------------------------------------------
template <class T>
struct std::tuple_size<quat::quaternion<T>> : std::integral_constant<std::size_t, 4> {};
template <std::size_t I, class T>
struct std::tuple_element<I, quat::quaternion<T>> {
    static_assert(I < 4, "quaternion index out of range");
    using type = T;
};

// ---- std::format support --------------------------------------------------------
// The format spec is that of T and is applied to each component:
//   std::format("{:.2f}", q)  ->  "(1.00,2.00,3.00,4.00)"
template <class T, class CharT>
struct std::formatter<quat::quaternion<T>, CharT> : std::formatter<T, CharT> {
    template <class FormatContext>
    auto format(const quat::quaternion<T>& q, FormatContext& ctx) const {
        auto out = ctx.out();
        const T c[4] = {q.w(), q.x(), q.y(), q.z()};
        *out++ = CharT('(');
        for (int n = 0; n < 4; ++n) {
            if (n) *out++ = CharT(',');
            ctx.advance_to(out);
            out = std::formatter<T, CharT>::format(c[n], ctx);
        }
        *out++ = CharT(')');
        return out;
    }
};

#endif  // QUATERNION_QUATERNION_HPP
