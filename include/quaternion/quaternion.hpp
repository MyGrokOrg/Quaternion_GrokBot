// quaternion.hpp - header-only C++23 quaternion numbers modeled on <complex>.
// SPDX-License-Identifier: MIT
//
// A quaternion q = a + b*i + c*j + d*k, with i^2 = j^2 = k^2 = ijk = -1.
// Multiplication is associative but NOT commutative.
#ifndef QUATERNION_QUATERNION_HPP
#define QUATERNION_QUATERNION_HPP

#include <algorithm>
#include <cmath>
#include <complex>
#include <concepts>
#include <cstddef>
#include <format>
#include <istream>
#include <limits>
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
// Mathematical constants as long double literals.  (std::numbers is avoided because
// some standard libraries, e.g. libc++ 19, give only double precision for long double.)
template <class T>
inline constexpr T pi_v = static_cast<T>(3.141592653589793238462643383279502884L);
template <class T>
inline constexpr T ln2_v = static_cast<T>(0.693147180559945309417232121458176568L);
template <class T>
inline constexpr T ln10_v = static_cast<T>(2.302585092994045684017991454684364208L);
template <class T>
inline constexpr T egamma_v = static_cast<T>(0.577215664901532860606512090082402431L);

// ---- C99 Annex G (G.5.1) recovery of infinities in multiplication and division ----
// Hamilton product r = p * q of raw components.
template <class T>
constexpr void hamilton(const T* p, const T* q, T* r) noexcept {
    r[0] = p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3];
    r[1] = p[0] * q[1] + p[1] * q[0] + p[2] * q[3] - p[3] * q[2];
    r[2] = p[0] * q[2] - p[1] * q[3] + p[2] * q[0] + p[3] * q[1];
    r[3] = p[0] * q[3] + p[1] * q[2] - p[2] * q[1] + p[3] * q[0];
}
template <class T>
inline bool all_nan(const T* r) noexcept {
    return std::isnan(r[0]) && std::isnan(r[1]) && std::isnan(r[2]) && std::isnan(r[3]);
}
template <class T>
inline bool any_inf(const T* p) noexcept {
    return std::isinf(p[0]) || std::isinf(p[1]) || std::isinf(p[2]) || std::isinf(p[3]);
}
template <class T>
inline bool all_finite(const T* p) noexcept {
    return std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2]) && std::isfinite(p[3]);
}
// "Box" an infinite operand: infinite components become +-1, the others +-0.
template <class T>
inline void box(T* p) noexcept {
    for (int i = 0; i < 4; ++i) p[i] = std::copysign(std::isinf(p[i]) ? T(1) : T(), p[i]);
}
template <class T>
inline void nan_to_zero(T* p) noexcept {
    for (int i = 0; i < 4; ++i)
        if (std::isnan(p[i])) p[i] = std::copysign(T(), p[i]);
}
// The _Cmultd recovery of C11 G.5.1, generalized to four components: applied when every
// component of the naive product is NaN.
template <class T>
inline void mul_recover(T* p, T* q, T* r) noexcept {
    if (!all_nan(r)) return;
    bool recalc = false;
    if (any_inf(p)) {
        box(p);
        nan_to_zero(q);
        recalc = true;
    }
    if (any_inf(q)) {
        box(q);
        nan_to_zero(p);
        recalc = true;
    }
    if (!recalc) {
        // overflow of a partial product
        bool ovf = false;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) ovf = ovf || std::isinf(p[i] * q[j]);
        if (ovf) {
            nan_to_zero(p);
            nan_to_zero(q);
            recalc = true;
        }
    }
    if (recalc) {
        hamilton(p, q, r);
        for (int i = 0; i < 4; ++i) r[i] *= std::numeric_limits<T>::infinity();
    }
}
// Division p / q = p conj(q) / |q|^2 with q scaled by a power of two (as _Cdivd), and the
// _Cdivd recovery when every component of the result is NaN.
template <class T>
inline void div_annexg(T* p, T* q, T* r) noexcept {
    const T m = std::max({std::abs(q[0]), std::abs(q[1]), std::abs(q[2]), std::abs(q[3])});
    const T logbw = std::logb(m);
    int k = 0;
    T qs[4] = {q[0], q[1], q[2], q[3]};
    if (std::isfinite(logbw)) {
        k = static_cast<int>(logbw);
        for (int i = 0; i < 4; ++i) qs[i] = std::scalbn(q[i], -k);
    }
    const T denom = qs[0] * qs[0] + qs[1] * qs[1] + qs[2] * qs[2] + qs[3] * qs[3];
    const T cq[4] = {qs[0], -qs[1], -qs[2], -qs[3]};
    T h[4];
    hamilton(p, cq, h);
    for (int i = 0; i < 4; ++i) r[i] = std::scalbn(h[i] / denom, -k);
    if (!all_nan(r)) return;
    const T inf = std::numeric_limits<T>::infinity();
    const bool p_nan = std::isnan(p[0]) && std::isnan(p[1]) && std::isnan(p[2]) && std::isnan(p[3]);
    if (denom == T() && !p_nan) {
        for (int i = 0; i < 4; ++i) r[i] = std::copysign(inf, q[0]) * p[i];
    } else if (any_inf(p) && all_finite(q)) {
        box(p);
        const T c2[4] = {q[0], -q[1], -q[2], -q[3]};
        hamilton(p, c2, h);
        for (int i = 0; i < 4; ++i) r[i] = inf * h[i];
    } else if (logbw == inf && all_finite(p)) {
        box(q);
        const T c2[4] = {q[0], -q[1], -q[2], -q[3]};
        hamilton(p, c2, h);
        for (int i = 0; i < 4; ++i) r[i] = T() * h[i];
    }
}
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
        T p[4] = {a_, b_, c_, d_};
        T r2[4] = {static_cast<T>(q.R_component_1()), static_cast<T>(q.R_component_2()),
                   static_cast<T>(q.R_component_3()), static_cast<T>(q.R_component_4())};
        T r[4];
        detail::hamilton(p, r2, r);
        if !consteval {
            detail::mul_recover(p, r2, r);  // Annex G infinity recovery (no-op if finite)
        }
        a_ = r[0]; b_ = r[1]; c_ = r[2]; d_ = r[3];
        return *this;
    }
    /// Right-division: *this = *this * q^-1.
    template <std::floating_point U>
    constexpr quaternion& operator/=(const quaternion<U>& q) noexcept {
        const quaternion r(q);
        if consteval {
            const T n = r.a_ * r.a_ + r.b_ * r.b_ + r.c_ * r.c_ + r.d_ * r.d_;
            *this *= quaternion(r.a_ / n, -r.b_ / n, -r.c_ / n, -r.d_ / n);
        } else {
            // p conj(q) / |q|^2 with q scaled by a power of two; Annex G recovery.
            T p[4] = {a_, b_, c_, d_}, d[4] = {r.a_, r.b_, r.c_, r.d_}, o[4];
            detail::div_annexg(p, d, o);
            a_ = o[0]; b_ = o[1]; c_ = o[2]; d_ = o[3];
        }
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
    // Nested two-argument hypot: libc++ 18's three-argument std::hypot squares its
    // arguments without scaling and overflows/underflows (e.g. for float beyond 2^64).
    return std::hypot(std::hypot(q.x(), q.y()), q.z());
}
// Rebuild s + u*v where u is the unit vector of q (i if q has no vector part).
template <class T>
[[nodiscard]] inline quaternion<T> along(const quaternion<T>& q, T vlen, T s, T v) noexcept {
    if (vlen == T()) return {s, v, T(), T()};
    const T f = v / vlen;
    // A zero component times an infinite f is a signed zero, not 0*inf = NaN.
    const auto m = [f](T x) { return x == T() ? std::copysign(T(), x) * (std::signbit(f) ? T(-1) : T(1)) : x * f; };
    return {s, m(q.x()), m(q.y()), m(q.z())};
}

/// Signed zero carrying the sign of the product i * u (copysign(0, u) if i is NaN).
template <class T>
[[nodiscard]] inline T zero_prod(T i, T u) noexcept {
    const bool neg = std::isnan(i) ? std::signbit(u) : (std::signbit(i) != std::signbit(u));
    return neg ? -T() : T();
}
template <class T>
[[nodiscard]] inline bool is_finite(const quaternion<T>& q) noexcept {
    return std::isfinite(q.w()) && std::isfinite(q.x()) && std::isfinite(q.y()) &&
           std::isfinite(q.z());
}
template <class T>
[[nodiscard]] inline bool has_nan(const quaternion<T>& q) noexcept {
    return std::isnan(q.w()) || std::isnan(q.x()) || std::isnan(q.y()) || std::isnan(q.z());
}

/// Special-value lift (C99/C11 Annex G mapped to quaternions by axis reduction).
/// f maps std::complex<T> -> std::complex<T> and is evaluated on the i-axis image z of q:
///  * vector part zero: z = a + i*b with b = x the *signed* i-component, so q = a + i*(+-0)
///    is exactly the std::complex argument (a, +-0); the j/k zeros of the result keep their
///    sign, flipped when the sign of Im f differs from the sign of b (the zeros of the vector
///    part behave as an infinitesimal vector parallel to i).
///  * otherwise z = a + i*|v| and the result is Re f + u * Im f with u the unit axis; for an
///    infinite v, u has components +-1/sqrt(m) on its m infinite components and signed zeros
///    elsewhere; for a NaN v (no infinity) the NaN components receive Im f.  Zero components
///    of u always yield signed zeros (sign of the product), never 0*inf = NaN.
template <class T, class F>
[[nodiscard]] inline quaternion<T> lift_sv(const quaternion<T>& q, F&& f) {
    const T vx = q.x(), vy = q.y(), vz = q.z();
    if (vx == T() && vy == T() && vz == T()) {
        const std::complex<T> r = f(std::complex<T>(q.w(), vx));
        const bool flip = !std::isnan(r.imag()) && std::signbit(r.imag()) != std::signbit(vx);
        const auto z0 = [flip](T z) { return flip ? -z : z; };
        return {r.real(), r.imag(), z0(vy), z0(vz)};
    }
    const T vl = vabs(q);
    const std::complex<T> r = f(std::complex<T>(q.w(), vl));
    const T I = r.imag();
    T u[3] = {vx, vy, vz};
    if (std::isinf(vl)) {
        int m = 0;
        for (T x : u) m += std::isinf(x);
        const T c = T(1) / std::sqrt(T(m));
        for (T& x : u) x = std::isinf(x) ? std::copysign(c, x) : std::isnan(x) ? x : std::copysign(T(), x);
    } else if (std::isnan(vl)) {
        for (T& x : u) x = std::isnan(x) ? T(1) : std::copysign(T(), x);
    } else {
        for (T& x : u) x /= vl;
    }
    T o[3];
    for (int j = 0; j < 3; ++j) o[j] = (u[j] == T()) ? zero_prod(I, u[j]) : I * u[j];
    return {r.real(), o[0], o[1], o[2]};
}
/// True when the special-value path must be used (non-finite input or zero vector part).
template <class T>
[[nodiscard]] inline bool needs_sv(const quaternion<T>& q) noexcept {
    return !is_finite(q) || (q.x() == T() && q.y() == T() && q.z() == T());
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
    const T c = std::cos(theta), s = std::sin(theta);
    // avoid inf * 0 = NaN for an infinite rho at theta = 0
    return {c == T() ? detail::zero_prod(rho, c) : rho * c,
            s == T() ? detail::zero_prod(rho, s) : rho * s, T(), T()};
}
/// rho * (cos(theta) + u sin(theta)) with u = axis/|axis| (the vector part of axis).
template <class T>
[[nodiscard]] inline quaternion<T> polar(const T& rho, const T& theta,
                                         const quaternion<T>& axis) noexcept {
    const T c = std::cos(theta), s = std::sin(theta);
    return detail::along(axis, detail::vabs(axis), c == T() ? detail::zero_prod(rho, c) : rho * c,
                         s == T() ? detail::zero_prod(rho, s) : rho * s);
}

// ---- transcendental functions ----------------------------------------------
// Each function has a finite-path evaluator (finite input, nonzero vector part) and an
// Annex G kernel ag::f(x, y) on std::complex that handles the special values (C11 G.6),
// reduced to the first quadrant by the symmetries of G.6 and evaluated on the finite
// path otherwise.  Public functions dispatch through detail::lift_sv for non-finite
// inputs and inputs with a zero vector part.
namespace detail::fin {
template <class T>
[[nodiscard]] inline quaternion<T> exp(const quaternion<T>& q) noexcept {
    const T v = vabs(q), e = std::exp(q.w());
    return along(q, v, e * std::cos(v), e * std::sin(v));
}
template <class T>
[[nodiscard]] inline quaternion<T> log(const quaternion<T>& q) noexcept {
    const T v = vabs(q);
    return along(q, v, std::log(quat::abs(q)), std::atan2(v, q.w()));
}
template <class T>
[[nodiscard]] inline quaternion<T> sqrt(const quaternion<T>& q) noexcept {
    // No cancellation: the larger of Re and |Im| comes from (r + |a|)/2, the other from
    // v / (2 * larger).  (r + |a|)/2 is formed as r/2 + |a|/2 to avoid overflow.
    const T v = vabs(q), r = quat::abs(q), a = q.w();
    const T big = std::sqrt(r / T(2) + std::abs(a) / T(2));
    if (big == T()) return along(q, v, T(), std::sqrt(r));
    const T small = (v / big) / T(2);
    return a >= T() ? along(q, v, big, small) : along(q, v, small, big);
}
template <class T>
[[nodiscard]] inline quaternion<T> sinh(const quaternion<T>& q) noexcept {
    const T v = vabs(q);
    return along(q, v, std::sinh(q.w()) * std::cos(v), std::cosh(q.w()) * std::sin(v));
}
template <class T>
[[nodiscard]] inline quaternion<T> cosh(const quaternion<T>& q) noexcept {
    const T v = vabs(q);
    return along(q, v, std::cosh(q.w()) * std::cos(v), std::sinh(q.w()) * std::sin(v));
}
// Kahan's complex tanh ("Branch cuts for complex elementary functions", 1987):
// tanh(x + i y) without overflow for any finite x, y.
template <class T>
[[nodiscard]] inline std::complex<T> kahan_tanh(T x, T y) noexcept {
    const T cut = T(std::numeric_limits<T>::digits + 4) * ln2_v<T> / T(2);
    if (std::abs(x) > cut) {
        const T e = std::exp(-std::abs(x));
        return {std::copysign(T(1), x), T(4) * std::sin(y) * std::cos(y) * e * e};
    }
    const T t = std::tan(y), b = T(1) + t * t, sh = std::sinh(x), rho = std::sqrt(T(1) + sh * sh);
    const T den = T(1) + b * sh * sh;
    return {b * rho * sh / den, t / den};
}
template <class T>
[[nodiscard]] inline quaternion<T> tanh(const quaternion<T>& q) noexcept {
    const T v = vabs(q);
    const std::complex<T> w = kahan_tanh(q.w(), v);
    return along(q, v, w.real(), w.imag());
}
template <class T>
[[nodiscard]] inline quaternion<T> sin(const quaternion<T>& q) noexcept {
    const T v = vabs(q);
    return along(q, v, std::sin(q.w()) * std::cosh(v), std::cos(q.w()) * std::sinh(v));
}
template <class T>
[[nodiscard]] inline quaternion<T> cos(const quaternion<T>& q) noexcept {
    const T v = vabs(q);
    return along(q, v, std::cos(q.w()) * std::cosh(v), -std::sin(q.w()) * std::sinh(v));
}
template <class T>
[[nodiscard]] inline quaternion<T> tan(const quaternion<T>& q) noexcept {
    // tan(a + u v) = -i tanh(-v + i a) on the plane of u
    const T v = vabs(q);
    const std::complex<T> w = kahan_tanh(-v, q.w());
    return along(q, v, w.imag(), -w.real());
}
}  // namespace detail::fin

namespace detail::ag {
template <class T>
using C = std::complex<T>;
template <class T>
[[nodiscard]] inline bool fin(T x) noexcept { return std::isfinite(x); }
// first-quadrant reduction helpers: y is +0, positive or NaN after these
template <class T>
[[nodiscard]] inline bool negy(T y) noexcept { return !std::isnan(y) && std::signbit(y); }
template <class T>
[[nodiscard]] inline C<T> cj(const C<T>& z) noexcept { return {z.real(), -z.imag()}; }
template <class T>
[[nodiscard]] inline C<T> nn() noexcept {
    return {std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN()};
}
// finite path at x + i y, y > 0
template <class T, class F>
[[nodiscard]] inline C<T> on_i(F f, T x, T y) noexcept {
    const quaternion<T> r = f(quaternion<T>(x, y, T(), T()));
    return {r.w(), r.x()};
}
// infinity * cis(y) for finite nonzero y (the exp/sinh/cosh entries)
template <class T>
[[nodiscard]] inline C<T> inf_cis(T y) noexcept {
    const T inf = std::numeric_limits<T>::infinity();
    return {inf * std::cos(y), inf * std::sin(y)};
}
template <class T>
[[nodiscard]] inline C<T> exp(T x, T y) noexcept {
    if (negy(y)) return cj(exp(x, -y));
    const T inf = std::numeric_limits<T>::infinity();
    if (fin(x)) {
        if (y == T()) return {std::exp(x), y};
        return fin(y) ? on_i([](const quaternion<T>& q) { return fin::exp(q); }, x, y) : nn<T>();
    }
    if (x == inf) return y == T() ? C<T>(inf, y) : fin(y) ? inf_cis(y) : C<T>(inf, y - y);
    if (x == -inf)
        return fin(y) ? C<T>(std::copysign(T(), std::cos(y)), std::copysign(T(), std::sin(y)))
                      : C<T>(T(), T());
    return y == T() ? C<T>(x, y) : nn<T>();
}
template <class T>
[[nodiscard]] inline C<T> log(T x, T y) noexcept {
    if (negy(y)) return cj(log(x, -y));
    const T inf = std::numeric_limits<T>::infinity(), pi = pi_v<T>;
    if (std::isinf(y)) return std::isnan(x) ? C<T>(inf, x) : C<T>(inf, std::atan2(y, x));
    if (std::isnan(y)) return std::isinf(x) ? C<T>(inf, y) : nn<T>();
    if (std::isinf(x)) return {inf, std::atan2(y, x)};
    if (std::isnan(x)) return nn<T>();
    if (y == T()) {
        if (x == T()) return {-inf, std::signbit(x) ? pi : T()};
        return x > T() ? C<T>(std::log(x), y) : C<T>(std::log(-x), pi);
    }
    return on_i([](const quaternion<T>& q) { return fin::log(q); }, x, y);
}
template <class T>
[[nodiscard]] inline C<T> log10(T x, T y) noexcept {
    const C<T> r = log(x, y);
    if (fin(x) && fin(y) && y != T()) {
        const quaternion<T> q = fin::log(quaternion<T>(x, std::abs(y), T(), T())) / ln10_v<T>;
        return {q.w(), std::copysign(q.x(), y)};
    }
    return {r.real() / ln10_v<T>, r.imag() / ln10_v<T>};
}
template <class T>
[[nodiscard]] inline C<T> sqrt(T x, T y) noexcept {
    if (negy(y)) return cj(sqrt(x, -y));
    const T inf = std::numeric_limits<T>::infinity();
    if (y == inf) return {inf, inf};
    if (x == -inf) return fin(y) ? C<T>(T(), inf) : C<T>(y, inf);
    if (x == inf) return fin(y) ? C<T>(inf, T()) : C<T>(inf, y);
    if (std::isnan(x) || std::isnan(y)) return nn<T>();
    if (y == T()) return x >= T() ? C<T>(x == T() ? T() : std::sqrt(x), y) : C<T>(T(), std::sqrt(-x));
    return on_i([](const quaternion<T>& q) { return fin::sqrt(q); }, x, y);
}
template <class T>
[[nodiscard]] inline C<T> sinh(T x, T y) noexcept {
    if (!std::isnan(x) && std::signbit(x)) return -sinh(-x, -y);
    if (negy(y)) return cj(sinh(x, -y));
    const T inf = std::numeric_limits<T>::infinity();
    if (fin(x)) {
        if (!fin(y)) return x == T() ? C<T>(x, y - y) : nn<T>();
        if (y == T()) return {std::sinh(x), y};
        return on_i([](const quaternion<T>& q) { return fin::sinh(q); }, x, y);
    }
    if (x == inf) return y == T() ? C<T>(x, y) : fin(y) ? inf_cis(y) : C<T>(inf, y - y);
    return y == T() ? C<T>(x, y) : nn<T>();
}
template <class T>
[[nodiscard]] inline C<T> cosh(T x, T y) noexcept {
    if (!std::isnan(x) && std::signbit(x)) return cosh(-x, -y);
    if (negy(y)) return cj(cosh(x, -y));
    const T inf = std::numeric_limits<T>::infinity();
    if (fin(x)) {
        if (!fin(y)) return x == T() ? C<T>(y - y, T()) : nn<T>();
        if (y == T()) return {std::cosh(x), y};
        return on_i([](const quaternion<T>& q) { return fin::cosh(q); }, x, y);
    }
    if (x == inf) return y == T() ? C<T>(x, y) : fin(y) ? inf_cis(y) : C<T>(inf, y - y);
    return y == T() ? C<T>(x, y) : nn<T>();
}
template <class T>
[[nodiscard]] inline C<T> tanh(T x, T y) noexcept {
    if (!std::isnan(x) && std::signbit(x)) return -tanh(-x, -y);
    if (negy(y)) return cj(tanh(x, -y));
    const T inf = std::numeric_limits<T>::infinity();
    if (fin(x)) {
        if (!fin(y)) return x == T() ? C<T>(x, y - y) : nn<T>();  // C23: ctanh(+0+i inf) = +0+iNaN
        if (y == T()) return {std::tanh(x), y};
        return on_i([](const quaternion<T>& q) { return fin::tanh(q); }, x, y);
    }
    if (x == inf) return {T(1), fin(y) ? std::copysign(T(), std::sin(T(2) * y)) : T()};
    return y == T() ? C<T>(x, y) : nn<T>();
}
// f(z) = -i g(i z): with g(-y + i x) = R + i I, f = I - i R (exact)
template <class T>
[[nodiscard]] inline C<T> mi(const C<T>& w) noexcept { return {w.imag(), -w.real()}; }
template <class T>
[[nodiscard]] inline C<T> sin(T x, T y) noexcept { return mi(sinh(-y, x)); }
template <class T>
[[nodiscard]] inline C<T> cos(T x, T y) noexcept { return cosh(-y, x); }
template <class T>
[[nodiscard]] inline C<T> tan(T x, T y) noexcept { return mi(tanh(-y, x)); }
}  // namespace detail::ag

#define QUAT_AG_FN(NAME)                                                                   \
    template <class T>                                                                     \
    [[nodiscard]] inline quaternion<T> NAME(const quaternion<T>& q) noexcept {             \
        if (detail::needs_sv(q))                                                           \
            return detail::lift_sv(                                                        \
                q, [](const std::complex<T>& z) { return detail::ag::NAME(z.real(), z.imag()); }); \
        return detail::fin::NAME(q);                                                       \
    }
/// exp(a + u v) = e^a (cos v + u sin v).
QUAT_AG_FN(exp)
/// Principal logarithm ln|q| + u*arg(q). For negative reals ln|a| + pi*i (-pi*i if the i
/// component is -0, as std::log).
QUAT_AG_FN(log)
/// Principal square root (non-negative real part).
QUAT_AG_FN(sqrt)
QUAT_AG_FN(sin)
QUAT_AG_FN(cos)
QUAT_AG_FN(tan)
QUAT_AG_FN(sinh)
QUAT_AG_FN(cosh)
QUAT_AG_FN(tanh)
#undef QUAT_AG_FN
template <class T>
[[nodiscard]] inline quaternion<T> log10(const quaternion<T>& q) noexcept {
    if (detail::needs_sv(q))
        return detail::lift_sv(
            q, [](const std::complex<T>& z) { return detail::ag::log10(z.real(), z.imag()); });
    return detail::fin::log(q) / detail::ln10_v<T>;
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
