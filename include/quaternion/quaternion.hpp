// quaternion.hpp - header-only C++23 quaternion numbers modeled on <complex>.
// SPDX-License-Identifier: MIT
//
// A quaternion q = a + b*i + c*j + d*k, with i^2 = j^2 = k^2 = ijk = -1.
// Multiplication is associative but NOT commutative.
#ifndef QUATERNION_QUATERNION_HPP
#define QUATERNION_QUATERNION_HPP

#include <algorithm>
#include <bit>
#include <cmath>
#include <complex>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <format>
#include <istream>
#include <iterator>
#include <limits>
#include <ostream>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

// QUATERNION_X87_LIBM: 1 = use the library's own long double sin, cos, tan, sinh, expm1, pow,
// tgamma and lgamma (quaternion/detail/x87_libm.hpp) instead of the C library's; needs the x87 80-bit long double
// (x86_64, GCC or Clang).  Off (0) by default for now; -DQUATERNION_X87_LIBM=1 forces it on, which
// CI does on Linux to test it.  See the README, "Configuration macros".
#ifndef QUATERNION_X87_LIBM
#define QUATERNION_X87_LIBM 0
#endif
#if QUATERNION_X87_LIBM
#include "detail/x87_libm.hpp"
#endif

// MSVC with /fp:contract or /fp:fast (both define _M_FP_CONTRACT) fuses a * b + c into an fma
// across statements and inlined calls, e.g. the product inside two_prod with the sum that
// follows it in fast_two_sum.  That breaks the error-free transformations (double-length
// arithmetic, exact squares) the library relies on, so the library's own functions are
// compiled with contraction off there.  Each header saves the includer's floating-point
// pragma state with float_control(push) and restores it with float_control(pop) at its end,
// so a user's #pragma fp_contract(off) (or on) before the include still holds after it.
// Explicit std::fma calls are unaffected.  GCC and Clang (including clang-cl) need nothing.
#if defined(_MSC_VER) && !defined(__clang__) && defined(_M_FP_CONTRACT)
#define QUATERNION_MSVC_CONTRACT_OFF 1
#pragma float_control(push)
#pragma fp_contract(off)
#endif

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

// Floating-point conversion U -> T is not narrowing.  For the standard types this is by
// conversion rank (so long double -> double is narrowing even where both are binary64, as
// for std::complex); otherwise by list-initialization.  (MSVC 19.5x evaluates the
// requires-expression form inside explicit() wrongly for float <- double, making the
// narrowing constructor unusable even explicitly.)
template <class T>
inline constexpr int fp_rank = std::is_same_v<T, float> ? 1 : std::is_same_v<T, double> ? 2
                             : std::is_same_v<T, long double> ? 3 : 0;
template <class U, class T>
consteval bool fp_widens() {
    if constexpr (fp_rank<U> != 0 && fp_rank<T> != 0) return fp_rank<U> <= fp_rank<T>;
    else return requires(U u) { T{u}; };
}

// ---- routing of the C library's transcendental functions (detail::lm) ----------
// The headers call these instead of std::sin, std::cos, ... (also for std::complex
// arguments).  Each forwards to the std:: function with the same arguments, so the overload
// chosen, the results and the code are those of the std:: call.  Where QUATERNION_X87_LIBM is
// on (see quaternion/detail/x87_libm.hpp), that header adds non-template overloads for x87
// long double, which overload resolution prefers.
namespace lm {
#define QUAT_LM_FWD(NAME)                                                                   \
    template <class... A>                                                                   \
    [[nodiscard]] constexpr auto NAME(const A&... a) noexcept(noexcept(std::NAME(a...)))    \
        -> decltype(std::NAME(a...)) {                                                      \
        return std::NAME(a...);                                                             \
    }
QUAT_LM_FWD(sin)
QUAT_LM_FWD(cos)
QUAT_LM_FWD(tan)
QUAT_LM_FWD(sinh)
QUAT_LM_FWD(cosh)
QUAT_LM_FWD(exp)
QUAT_LM_FWD(expm1)
QUAT_LM_FWD(pow)
QUAT_LM_FWD(tgamma)
QUAT_LM_FWD(lgamma)
QUAT_LM_FWD(polar)
QUAT_LM_FWD(sqrt)
#undef QUAT_LM_FWD
}  // namespace lm

// ---- constexpr floating-point helpers -------------------------------------------
// <cmath> is not constexpr on Clang or libc++ before C++26 (P1383), so the
// transcendental-free operations use these in constant evaluation.  Every runtime branch
// calls the same std:: function as before, so runtime results and code are unchanged.
// Special values in constant evaluation: operations that would produce a NaN, or (on GCC)
// overflow, are not constant expressions, so only the finite paths need to be exact here.
namespace cx {
template <class T>
constexpr bool isnan(T x) noexcept { return x != x; }
template <class T>
constexpr bool isinf(T x) noexcept {
    return x == std::numeric_limits<T>::infinity() || x == -std::numeric_limits<T>::infinity();
}
template <class T>
constexpr bool isfinite(T x) noexcept { return !isnan(x) && !isinf(x); }

// Sign bit (also of zeros and NaNs) through std::bit_cast of the representation.
template <class T>
constexpr bool signbit(T x) noexcept {
    using L = std::numeric_limits<T>;
    if constexpr (sizeof(T) == 4 && L::digits == 24) {
        return (std::bit_cast<std::uint32_t>(x) >> 31) != 0;
    } else if constexpr (sizeof(T) == 8 && L::digits == 53) {
        return (std::bit_cast<std::uint64_t>(x) >> 63) != 0;
    } else if constexpr (L::digits == 64 && sizeof(T) == 16) {  // x87 80-bit in 16 bytes
        struct rep { std::uint64_t m; std::uint16_t se; };
        static_assert(std::endian::native == std::endian::little);
        return (std::bit_cast<rep>(x).se >> 15) != 0;
    } else if constexpr (L::digits == 113 && sizeof(T) == 16) {  // IEEE binary128
        struct rep { std::uint64_t w[2]; };
        const rep r = std::bit_cast<rep>(x);
        return (r.w[std::endian::native == std::endian::little ? 1 : 0] >> 63) != 0;
    } else {
        // Other formats: exact for everything but negative NaNs.
        return x < T() || (x == T() && T(1) / x < T());
    }
}
template <class T>
constexpr T fabs(T x) noexcept { return signbit(x) ? -x : x; }
template <class T>
constexpr T copysign(T x, T s) noexcept { return signbit(s) ? -fabs(x) : fabs(x); }

// Representation access for the formats with a known layout (O(1), so constant
// evaluation stays far below the compilers' step limits).
template <class T>
inline constexpr int fmt_kind = (sizeof(T) == 4 && std::numeric_limits<T>::digits == 24) ? 1
                              : (sizeof(T) == 8 && std::numeric_limits<T>::digits == 53) ? 2
                              : (sizeof(T) == 16 && std::numeric_limits<T>::digits == 64 &&
                                 std::endian::native == std::endian::little) ? 3
                              : 0;
struct x87_rep { std::uint64_t m; std::uint16_t se; std::uint16_t pad1; std::uint32_t pad2; };  // for writing
struct x87_view { std::uint64_t m; std::uint16_t se; };  // for reading (the 6 padding bytes stay padding)

// 2^e for any e from the smallest subnormal exponent to max_exponent - 1 (exact).
template <class T>
constexpr T pow2(int e) noexcept {
    using L = std::numeric_limits<T>;
    constexpr int dmin = L::min_exponent - L::digits;  // exponent of denorm_min
    if constexpr (fmt_kind<T> == 1) {
        return std::bit_cast<T>(e >= L::min_exponent - 1 ? static_cast<std::uint32_t>(e + 127) << 23
                                                         : std::uint32_t(1) << (e - dmin));
    } else if constexpr (fmt_kind<T> == 2) {
        return std::bit_cast<T>(e >= L::min_exponent - 1 ? static_cast<std::uint64_t>(e + 1023) << 52
                                                         : std::uint64_t(1) << (e - dmin));
    } else if constexpr (fmt_kind<T> == 3) {
        return e >= L::min_exponent - 1
                   ? std::bit_cast<T>(x87_rep{std::uint64_t(1) << 63, static_cast<std::uint16_t>(e + 16383), 0, 0})
                   : std::bit_cast<T>(x87_rep{std::uint64_t(1) << (e - dmin), 0, 0, 0});
    } else {
        if (e < L::min_exponent - 1) return pow2<T>(L::min_exponent - 1) * pow2<T>(e - (L::min_exponent - 1));
        T r(1), b = e < 0 ? T(0.5) : T(2);
        for (unsigned k = static_cast<unsigned>(e < 0 ? -e : e); k; k >>= 1) {
            if (k & 1u) r *= b;
            if (k > 1u) b *= b;  // (squaring past the last bit could overflow)
        }
        return r;
    }
}
// Exponent of a finite nonzero x: x = m 2^e with 1 <= |m| < 2 (as std::ilogb).
template <class T>
constexpr int ilogb(T x) noexcept {
    using L = std::numeric_limits<T>;
    constexpr int dmin = L::min_exponent - L::digits;
    if constexpr (fmt_kind<T> == 1 || fmt_kind<T> == 2) {
        using U = std::conditional_t<fmt_kind<T> == 1, std::uint32_t, std::uint64_t>;
        constexpr int mbits = L::digits - 1;
        const U u = std::bit_cast<U>(x);
        const int be = static_cast<int>((u >> mbits) & ((U(1) << (sizeof(U) * 8 - 1 - mbits)) - 1));
        if (be) return be - (L::max_exponent - 1);
        return std::bit_width(u & ((U(1) << mbits) - 1)) - 1 + dmin;  // subnormal
    } else if constexpr (fmt_kind<T> == 3) {
        const x87_view r = std::bit_cast<x87_view>(x);
        const int be = r.se & 0x7fff;
        if (be) return be - 16383;
        return std::bit_width(r.m) - 1 + dmin;
    } else {
        T a = fabs(x);
        int e = 0;
        if (a < L::min()) { a *= pow2<T>(L::digits); e -= L::digits; }
        for (int k = 1 << 14; k; k >>= 1)
            if (k < L::max_exponent && a >= pow2<T>(k)) { a *= pow2<T>(-k); e += k; }
        for (int k = 1 << 14; k; k >>= 1)
            if (k <= -L::min_exponent && a < pow2<T>(1 - k)) { a *= pow2<T>(k); e -= k; }
        while (a >= T(2)) { a /= T(2); ++e; }
        while (a < T(1)) { a *= T(2); --e; }
        return e;
    }
}
template <class T>
constexpr T logb(T x) noexcept {
    if (isnan(x)) return x;
    if (isinf(x)) return std::numeric_limits<T>::infinity();
    if (x == T()) return -std::numeric_limits<T>::infinity();
    return static_cast<T>(ilogb(x));
}
// x 2^n with a single rounding (as std::scalbn / std::ldexp in round-to-nearest).
template <class T>
constexpr T scalbn(T x, int n) noexcept {
    using L = std::numeric_limits<T>;
    if (x == T() || !isfinite(x) || n == 0) return x;
    const int e = ilogb(x);
    const int half = e / 2;
    T m = x * pow2<T>(-half) * pow2<T>(half - e);  // exact: 1 <= |m| < 2
    const long t = static_cast<long>(e) + n;           // target exponent
    if (t > L::max_exponent - 1) return copysign(L::infinity(), x);
    const int dmin = L::min_exponent - L::digits;      // exponent of denorm_min
    if (t < dmin - 1) return copysign(T(), x);
    if (t == dmin - 1) { m /= T(2); return m * pow2<T>(dmin); }
    return m * pow2<T>(static_cast<int>(t));          // exact if normal, one rounding if not
}
template <class T>
constexpr T ldexp(T x, int n) noexcept { return scalbn(x, n); }
// x = m 2^e with 0.5 <= |m| < 1 (as std::frexp).
template <class T>
constexpr T frexp(T x, int* e) noexcept {
    if (x == T() || !isfinite(x)) { *e = 0; return x; }
    const int k = ilogb(x) + 1;
    *e = k;
    return scalbn(x, -k);
}

// Error-free transformations (exact in constant evaluation, which never contracts).
template <class T>
constexpr void two_sum(T a, T b, T& s, T& e) noexcept {
    s = a + b;
    const T bb = s - a;
    e = (a - (s - bb)) + (b - bb);
}
template <class T>
constexpr void two_prod(T a, T b, T& p, T& e) noexcept {
    constexpr T split = pow2<T>((std::numeric_limits<T>::digits + 1) / 2) + T(1);
    p = a * b;
    const T ca = split * a, ah = ca - (ca - a), al = a - ah;
    const T cb = split * b, bh = cb - (cb - b), bl = b - bh;
    e = ((ah * bh - p) + ah * bl + al * bh) + al * bl;
}
// Correctly rounded square root of m given as the exact sum hi + lo, 1 <= hi < 4.
template <class T>
constexpr T sqrt_scaled(T hi, T lo) noexcept {
    // Newton from a crude start, then fix the last ulp with an exact midpoint test.
    T r = hi < T(2) ? T(1.25) : T(1.75);
    for (int i = 0; i < 8; ++i) r = (r + hi / r) / T(2);
    const T ulp = pow2<T>(ilogb(r) - std::numeric_limits<T>::digits + 1);
    for (int pass = 0; pass < 3; ++pass) {
        // sign of m - (r + s ulp/2)^2 for s = +1 and -1, with every term exact.
        const auto side = [&](T c) {  // sign of m - c^2 - c*ulp - ulp^2/4
            T p, pe, s1, e1, s2, e2, s3, e3, s4, e4;
            two_prod(c, c, p, pe);
            two_sum(hi, -p, s1, e1);
            two_sum(s1, lo - pe, s2, e2);
            two_sum(s2, -c * ulp, s3, e3);
            two_sum(s3, -(ulp * ulp) / T(4), s4, e4);
            const T tail = ((e1 + e2) + e3) + e4;
            const T v = s4 + tail;
            return v > T() ? 1 : v < T() ? -1 : 0;
        };
        if (side(r) > 0) { r += ulp; continue; }            // m > (r + ulp/2)^2
        if (side(r - ulp) < 0) { r -= ulp; continue; }      // m < (r - ulp/2)^2
        break;
    }
    return r;
}
template <class T>
constexpr T sqrt(T x) noexcept {
    if (isnan(x) || x == T() || x == std::numeric_limits<T>::infinity()) return x;
    if (x < T()) return std::numeric_limits<T>::quiet_NaN();
    int e = ilogb(x);
    if (e & 1) --e;  // even exponent: x = m 2^e, 1 <= m < 4
    const T m = scalbn(x, -e);
    return scalbn(sqrt_scaled(m, T()), e / 2);
}
// hypot(x, y) = sqrt(x^2 + y^2) without overflow or underflow, correctly rounded.
template <class T>
constexpr T hypot(T x, T y) noexcept {
    if (isinf(x) || isinf(y)) return std::numeric_limits<T>::infinity();
    if (isnan(x) || isnan(y)) return std::numeric_limits<T>::quiet_NaN();
    x = fabs(x); y = fabs(y);
    if (x < y) { const T t = x; x = y; y = t; }
    if (x == T()) return T();
    int e = ilogb(x);
    if (e & 1) --e;
    const T a = scalbn(x, -e), b = scalbn(y, -e);  // 1 <= a < 4; b may round if tiny
    T p, pe, q, qe, s, se;
    two_prod(a, a, p, pe);
    two_prod(b, b, q, qe);
    two_sum(p, q, s, se);
    T lo = se + (pe + qe);
    T hi = s;
    two_sum(hi, lo, hi, lo);
    // hi + lo in [1, 32): renormalize to [1, 4) by an exact even power of two.
    int f = 0;
    while (hi >= T(4)) { hi /= T(4); lo /= T(4); f += 1; }
    return scalbn(sqrt_scaled(hi, lo), e + f);
}
}  // namespace cx

// Wrappers: the constexpr helper in constant evaluation, the std:: call at run time.
template <class T> constexpr bool c_isnan(T x) noexcept { if consteval { return cx::isnan(x); } else { return std::isnan(x); } }
template <class T> constexpr bool c_isinf(T x) noexcept { if consteval { return cx::isinf(x); } else { return std::isinf(x); } }
template <class T> constexpr bool c_isfinite(T x) noexcept { if consteval { return cx::isfinite(x); } else { return std::isfinite(x); } }
template <class T> constexpr bool c_signbit(T x) noexcept { if consteval { return cx::signbit(x); } else { return std::signbit(x); } }
template <class T> constexpr T c_fabs(T x) noexcept { if consteval { return cx::fabs(x); } else { return std::abs(x); } }
template <class T> constexpr T c_copysign(T x, T s) noexcept { if consteval { return cx::copysign(x, s); } else { return std::copysign(x, s); } }
template <class T> constexpr T c_logb(T x) noexcept { if consteval { return cx::logb(x); } else { return std::logb(x); } }
// Apple's x86_64 libm rounds subnormal results of scalbn/scalbnl wrongly (236 of 400,000
// double and 76 of 400,000 long double inputs with subnormal results, macOS 15 x86_64), which
// changed the bits of quaternion division there; it gets the exact helper at run time too.
#if defined(__APPLE__) && (defined(__x86_64__) || defined(__i386__))
template <class T> constexpr T c_scalbn(T x, int n) noexcept { return cx::scalbn(x, n); }
#else
template <class T> constexpr T c_scalbn(T x, int n) noexcept { if consteval { return cx::scalbn(x, n); } else { return std::scalbn(x, n); } }
#endif
template <class T> constexpr T c_hypot(T x, T y) noexcept { if consteval { return cx::hypot(x, y); } else { return std::hypot(x, y); } }

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
constexpr bool all_nan(const T* r) noexcept {
    return c_isnan(r[0]) && c_isnan(r[1]) && c_isnan(r[2]) && c_isnan(r[3]);
}
template <class T>
constexpr bool any_inf(const T* p) noexcept {
    return c_isinf(p[0]) || c_isinf(p[1]) || c_isinf(p[2]) || c_isinf(p[3]);
}
template <class T>
constexpr bool all_finite(const T* p) noexcept {
    return c_isfinite(p[0]) && c_isfinite(p[1]) && c_isfinite(p[2]) && c_isfinite(p[3]);
}
// "Box" an infinite operand: infinite components become +-1, the others +-0.
template <class T>
constexpr void box(T* p) noexcept {
    for (int i = 0; i < 4; ++i) p[i] = c_copysign(c_isinf(p[i]) ? T(1) : T(), p[i]);
}
template <class T>
constexpr void nan_to_zero(T* p) noexcept {
    for (int i = 0; i < 4; ++i)
        if (c_isnan(p[i])) p[i] = c_copysign(T(), p[i]);
}
// The _Cmultd recovery of C11 G.5.1, generalized to four components: applied when every
// component of the naive product is NaN.
template <class T>
constexpr void mul_recover(T* p, T* q, T* r) noexcept {
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
            for (int j = 0; j < 4; ++j) ovf = ovf || c_isinf(p[i] * q[j]);
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
constexpr void div_annexg(T* p, T* q, T* r) noexcept {
    const T m = std::max({c_fabs(q[0]), c_fabs(q[1]), c_fabs(q[2]), c_fabs(q[3])});
    const T logbw = c_logb(m);
    int k = 0;
    T qs[4] = {q[0], q[1], q[2], q[3]};
    if (c_isfinite(logbw)) {
        k = static_cast<int>(logbw);
        for (int i = 0; i < 4; ++i) qs[i] = c_scalbn(q[i], -k);
    }
    const T denom = qs[0] * qs[0] + qs[1] * qs[1] + qs[2] * qs[2] + qs[3] * qs[3];
    const T cq[4] = {qs[0], -qs[1], -qs[2], -qs[3]};
    T h[4];
    hamilton(p, cq, h);
    for (int i = 0; i < 4; ++i) r[i] = c_scalbn(h[i] / denom, -k);
    if (!all_nan(r)) return;
    const T inf = std::numeric_limits<T>::infinity();
    const bool p_nan = c_isnan(p[0]) && c_isnan(p[1]) && c_isnan(p[2]) && c_isnan(p[3]);
    if (denom == T() && !p_nan) {
        for (int i = 0; i < 4; ++i) r[i] = c_copysign(inf, q[0]) * p[i];
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
    constexpr explicit(!detail::fp_widens<U, T>())
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
        // p conj(q) / |q|^2 with q scaled by a power of two; Annex G recovery.  The same
        // steps run in constant evaluation (with the detail::cx helpers), so a constexpr
        // quotient equals the runtime one wherever the runtime code is not contracted.
        T p[4] = {a_, b_, c_, d_}, d[4] = {r.a_, r.b_, r.c_, r.d_}, o[4];
        detail::div_annexg(p, d, o);
        a_ = o[0]; b_ = o[1]; c_ = o[2]; d_ = o[3];
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

    // ---- mixed precision: a floating-point scalar or complex of another type U.
    // The result is quaternion<common_type_t<T, U>>, the floating-point promotion that
    // [cmplx.over] uses (so quaternion<float> * 1.0 is a quaternion<double>, not a narrowed
    // quaternion<float>).  Integers keep converting to T through the operators above.
    // (std::complex itself has no mixed operators: complex<float>{} * 1.0 does not compile.)
#define QUAT_MIXOP(OP)                                                                       \
    template <std::floating_point U>                                                         \
        requires(!std::same_as<U, T>)                                                        \
    [[nodiscard]] friend constexpr quaternion<std::common_type_t<T, U>> operator OP(         \
        const quaternion& l, const U& r) noexcept {                                          \
        using R = std::common_type_t<T, U>;                                                  \
        return quaternion<R>(l) OP static_cast<R>(r);                                        \
    }                                                                                        \
    template <std::floating_point U>                                                         \
        requires(!std::same_as<U, T>)                                                        \
    [[nodiscard]] friend constexpr quaternion<std::common_type_t<T, U>> operator OP(         \
        const U& l, const quaternion& r) noexcept {                                          \
        using R = std::common_type_t<T, U>;                                                  \
        return static_cast<R>(l) OP quaternion<R>(r);                                        \
    }                                                                                        \
    template <std::floating_point U>                                                         \
        requires(!std::same_as<U, T>)                                                        \
    [[nodiscard]] friend constexpr quaternion<std::common_type_t<T, U>> operator OP(         \
        const quaternion& l, const std::complex<U>& r) noexcept {                            \
        using R = std::common_type_t<T, U>;                                                  \
        return quaternion<R>(l) OP std::complex<R>(static_cast<R>(r.real()), static_cast<R>(r.imag())); \
    }                                                                                        \
    template <std::floating_point U>                                                         \
        requires(!std::same_as<U, T>)                                                        \
    [[nodiscard]] friend constexpr quaternion<std::common_type_t<T, U>> operator OP(         \
        const std::complex<U>& l, const quaternion& r) noexcept {                            \
        using R = std::common_type_t<T, U>;                                                  \
        return std::complex<R>(static_cast<R>(l.real()), static_cast<R>(l.imag())) OP quaternion<R>(r); \
    }
    QUAT_MIXOP(+)
    QUAT_MIXOP(-)
    QUAT_MIXOP(*)
    QUAT_MIXOP(/)
#undef QUAT_MIXOP

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
    template <std::floating_point U>
        requires(!std::same_as<U, T>)
    [[nodiscard]] friend constexpr bool operator==(const quaternion& l, const U& r) noexcept {
        using R = std::common_type_t<T, U>;
        return quaternion<R>(l) == static_cast<R>(r);
    }
    template <std::floating_point U>
        requires(!std::same_as<U, T>)
    [[nodiscard]] friend constexpr bool operator==(const quaternion& l, const std::complex<U>& r) noexcept {
        using R = std::common_type_t<T, U>;
        return quaternion<R>(l) == std::complex<R>(static_cast<R>(r.real()), static_cast<R>(r.imag()));
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
// Windows and macOS: the C library's hypot is not correctly rounded (Universal CRT: 1 ulp off
// in 32576 of 200000 double inputs in test_constexpr; macOS: 31923 on arm64, 27884 on x86_64,
// and long double hypotl on x86_64 in 33178, some by more than 1 ulp), and every function of a
// quaternion with |Im| = |v| large inherits |v| times that error (e.g. sinh at |v| = 109:
// 63.8 eps on Windows).  There ||v|| is computed by vnorm3 below; glibc keeps std::hypot.
#if defined(_WIN32) || defined(__APPLE__)
#define QUATERNION_OWN_VNORM 1
#endif
#if defined(QUATERNION_OWN_VNORM)
// sqrt(x^2 + y^2 + z^2) to within about 0.5 ulp, exact when the true value is representable.
template <class T>
[[nodiscard]] inline T vnorm3(T x, T y, T z) noexcept {
    using L = std::numeric_limits<T>;
    x = std::abs(x); y = std::abs(y); z = std::abs(z);
    const T m = (std::max)((std::max)(x, y), z);
    if (!(m > T()) || !(m <= L::max()) || std::isnan(x) || std::isnan(y) || std::isnan(z))
        return std::hypot(std::hypot(x, y), z);  // zero, infinities and NaNs as before
    if constexpr (L::digits < 53 && std::numeric_limits<double>::digits == 53 && L::max_exponent <= 128) {
        // float: squares exact in double, sum rounded once, then sqrt and the final rounding.
        const double dx = x, dy = y, dz = z;
        return static_cast<T>(detail::lm::sqrt(dx * dx + dy * dy + dz * dz));
    } else {
        // Scale by 2^-e so that m is in [1, 2): no overflow, and squares of components below
        // m 2^-(digits+2) only matter beyond the last bit.
        const int e = std::ilogb(m);
        x = c_scalbn(x, -e); y = c_scalbn(y, -e); z = c_scalbn(z, -e);
        // Exact squares (Dekker), their sum in double-length, one Newton correction.
        const auto sq = [](T a, T& p, T& err) {
            constexpr T split = T((std::uint64_t(1) << ((L::digits + 1) / 2)) + 1);
            p = a * a;
            const T c = split * a, ah = c - (c - a), al = a - ah;
            err = ((ah * ah - p) + T(2) * ah * al) + al * al;
        };
        const auto two_sum = [](T a, T b, T& s2, T& err) {
            s2 = a + b;
            const T bb = s2 - a;
            err = (a - (s2 - bb)) + (b - bb);
        };
        T px, ex, py, ey, pz, ez, s1, e1, s2, e2;
        sq(x, px, ex); sq(y, py, ey); sq(z, pz, ez);
        two_sum(px, py, s1, e1);
        two_sum(s1, pz, s2, e2);
        const T lo = ((ex + ey) + ez) + (e1 + e2);
        const T h = detail::lm::sqrt(s2);
        // residual s - h^2, exactly: h^2 = hp + he
        T hp, he;
        sq(h, hp, he);
        const T r = ((s2 - hp) - he + lo) / (T(2) * h);
        return c_scalbn(h + r, e);
    }
}
#endif
template <class T>
[[nodiscard]] inline T vabs(const quaternion<T>& q) noexcept {
#if defined(QUATERNION_OWN_VNORM)
    return vnorm3(q.x(), q.y(), q.z());
#else
    // Nested two-argument hypot: libc++ 18's three-argument std::hypot squares its
    // arguments without scaling and overflows/underflows (e.g. for float beyond 2^64).
    return std::hypot(std::hypot(q.x(), q.y()), q.z());
#endif
}
// a * b = p + e exactly.  Hardware fma where the build has it (there the compiler may also
// contract the Dekker split, which would break it); Dekker's split otherwise and for
// formats wider than double (x87 long double has no fma instruction to contract into).
template <class T>
inline void exact_mul(T a, T b, T& p, T& e) noexcept {
    using L = std::numeric_limits<T>;
#if defined(FP_FAST_FMA) || defined(__FMA__) || defined(__AVX2__) || defined(__aarch64__) || defined(_M_ARM64)
    if constexpr (L::digits <= 53) {
        p = a * b;
        e = std::fma(a, b, -p);
        return;
    } else
#endif
    {
        constexpr T split = T((std::uint64_t(1) << ((L::digits + 1) / 2)) + 1);
        p = a * b;
        const T ca = split * a, ah = ca - (ca - a), al = a - ah;
        const T cb = split * b, bh = cb - (cb - b), bl = b - bh;
        e = ((ah * bh - p) + ah * bl + al * bh) + al * bl;
    }
}
// log|q| for 0.49 <= |q|^2 <= 2.02 (no overflow or harmful underflow in the squares), with
// relative accuracy: (log1p s)/2 with s = a^2 + b^2 + c^2 + d^2 - 1 rounded once from its
// exact value.  The eight parts of the exact squares and -1 are summed by repeated
// error-free cascades (Ogita-Rump-Oishi VecSum) until the leading term carries the sum.
template <class T>
[[nodiscard]] inline T log_abs_near1(const quaternion<T>& q) noexcept {
    T t[9];
    t[0] = T(-1);
    exact_mul(q.w(), q.w(), t[1], t[5]);
    exact_mul(q.x(), q.x(), t[2], t[6]);
    exact_mul(q.y(), q.y(), t[3], t[7]);
    exact_mul(q.z(), q.z(), t[4], t[8]);
    constexpr T eps = std::numeric_limits<T>::epsilon();
    for (int pass = 0; pass < 12; ++pass) {
        // cascade: t[8] becomes the rounded sum, t[0..7] the exact rounding errors
        for (int i = 1; i < 9; ++i) {
            const T a = t[i - 1], b = t[i], s = a + b, bb = s - a;
            t[i] = s;
            t[i - 1] = (a - (s - bb)) + (b - bb);
        }
        T tail = T();
        for (int i = 0; i < 8; ++i) tail += std::abs(t[i]);
        // the plain sum of the remaining tail then adds at most about 8 eps^2 |s|
        if (tail <= eps * std::abs(t[8])) break;
    }
    T lo = T();
    for (int i = 0; i < 8; ++i) lo += t[i];
    return std::log1p(t[8] + lo) / T(2);
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
        const T c = T(1) / detail::lm::sqrt(T(m));
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
[[nodiscard]] constexpr T abs(const quaternion<T>& q) noexcept {
    // Run time: std::hypot.  Constant evaluation: a correctly rounded constexpr hypot.
    using detail::c_hypot;
    return c_hypot(c_hypot(q.w(), q.x()), c_hypot(q.y(), q.z()));
}
/// Angle in [0, pi] between q and the positive real axis: atan2(|v|, a).
template <class T>
[[nodiscard]] inline T arg(const quaternion<T>& q) noexcept {
    return std::atan2(detail::vabs(q), q.w());
}
/// Unit quaternion q/|q| (q itself if zero).
template <class T>
[[nodiscard]] constexpr quaternion<T> normalize(const quaternion<T>& q) noexcept {
    const T r = abs(q);
    return r == T() ? q : q / r;
}
/// Projection onto the Riemann-sphere analogue (as std::proj).
template <class T>
[[nodiscard]] constexpr quaternion<T> proj(const quaternion<T>& q) noexcept {
    using detail::c_isinf, detail::c_copysign;
    if (c_isinf(q.w()) || c_isinf(q.x()) || c_isinf(q.y()) || c_isinf(q.z()))
        return {std::numeric_limits<T>::infinity(), c_copysign(T(), q.x()),
                c_copysign(T(), q.y()), c_copysign(T(), q.z())};
    return q;
}

/// rho * (cos(theta) + i sin(theta)) -- same as std::polar embedded on the i axis.
template <class T>
[[nodiscard]] inline quaternion<T> polar(const T& rho, const T& theta = T()) noexcept {
    const T c = detail::lm::cos(theta), s = detail::lm::sin(theta);
    // avoid inf * 0 = NaN for an infinite rho at theta = 0
    return {c == T() ? detail::zero_prod(rho, c) : rho * c,
            s == T() ? detail::zero_prod(rho, s) : rho * s, T(), T()};
}
/// rho * (cos(theta) + u sin(theta)) with u = axis/|axis| (the vector part of axis).
template <class T>
[[nodiscard]] inline quaternion<T> polar(const T& rho, const T& theta,
                                         const quaternion<T>& axis) noexcept {
    const T c = detail::lm::cos(theta), s = detail::lm::sin(theta);
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
    const T v = vabs(q), e = detail::lm::exp(q.w());
    return along(q, v, e * detail::lm::cos(v), e * detail::lm::sin(v));
}
template <class T>
[[nodiscard]] inline quaternion<T> log(const quaternion<T>& q) noexcept {
    const T v = vabs(q);
    // Next to |q| = 1, log(abs(q)) cancels (abs(q) carries an absolute error of ulp(1)/2);
    // there log|q| comes from |q|^2 - 1 in double length.  The plain sum of squares only
    // selects the path (overflow gives inf and underflow 0, both outside the range).
    const T n = q.w() * q.w() + q.x() * q.x() + q.y() * q.y() + q.z() * q.z();
    const T lr = (n >= T(0.49) && n <= T(2.02)) ? log_abs_near1(q) : std::log(quat::abs(q));
    return along(q, v, lr, std::atan2(v, q.w()));
}
template <class T>
[[nodiscard]] inline quaternion<T> sqrt(const quaternion<T>& q) noexcept {
    // No cancellation: the larger of Re and |Im| comes from (r + |a|)/2, the other from
    // v / (2 * larger).  (r + |a|)/2 is formed as r/2 + |a|/2 to avoid overflow.
    const T v = vabs(q), r = quat::abs(q), a = q.w();
    const T big = detail::lm::sqrt(r / T(2) + std::abs(a) / T(2));
    if (big == T()) return along(q, v, T(), detail::lm::sqrt(r));
    const T small = (v / big) / T(2);
    return a >= T() ? along(q, v, big, small) : along(q, v, small, big);
}
template <class T>
[[nodiscard]] inline quaternion<T> sinh(const quaternion<T>& q) noexcept {
    const T v = vabs(q);
    return along(q, v, detail::lm::sinh(q.w()) * detail::lm::cos(v), detail::lm::cosh(q.w()) * detail::lm::sin(v));
}
template <class T>
[[nodiscard]] inline quaternion<T> cosh(const quaternion<T>& q) noexcept {
    const T v = vabs(q);
    return along(q, v, detail::lm::cosh(q.w()) * detail::lm::cos(v), detail::lm::sinh(q.w()) * detail::lm::sin(v));
}
// Kahan's complex tanh ("Branch cuts for complex elementary functions", 1987):
// tanh(x + i y) without overflow for any finite x, y.
template <class T>
[[nodiscard]] inline std::complex<T> kahan_tanh(T x, T y) noexcept {
    const T cut = T(std::numeric_limits<T>::digits + 4) * ln2_v<T> / T(2);
    if (std::abs(x) > cut) {
        const T e = detail::lm::exp(-std::abs(x));
        return {std::copysign(T(1), x), T(4) * detail::lm::sin(y) * detail::lm::cos(y) * e * e};
    }
    const T t = detail::lm::tan(y), b = T(1) + t * t, sh = detail::lm::sinh(x), rho = detail::lm::sqrt(T(1) + sh * sh);
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
    return along(q, v, detail::lm::sin(q.w()) * detail::lm::cosh(v), detail::lm::cos(q.w()) * detail::lm::sinh(v));
}
template <class T>
[[nodiscard]] inline quaternion<T> cos(const quaternion<T>& q) noexcept {
    const T v = vabs(q);
    return along(q, v, detail::lm::cos(q.w()) * detail::lm::cosh(v), -detail::lm::sin(q.w()) * detail::lm::sinh(v));
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
    return {inf * detail::lm::cos(y), inf * detail::lm::sin(y)};
}
template <class T>
[[nodiscard]] inline C<T> exp(T x, T y) noexcept {
    if (negy(y)) return cj(exp(x, -y));
    const T inf = std::numeric_limits<T>::infinity();
    if (fin(x)) {
        if (y == T()) return {detail::lm::exp(x), y};
        return fin(y) ? on_i([](const quaternion<T>& q) { return fin::exp(q); }, x, y) : nn<T>();
    }
    if (x == inf) return y == T() ? C<T>(inf, y) : fin(y) ? inf_cis(y) : C<T>(inf, y - y);
    if (x == -inf)
        return fin(y) ? C<T>(std::copysign(T(), detail::lm::cos(y)), std::copysign(T(), detail::lm::sin(y)))
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
    if (y == T()) return x >= T() ? C<T>(x == T() ? T() : detail::lm::sqrt(x), y) : C<T>(T(), detail::lm::sqrt(-x));
    return on_i([](const quaternion<T>& q) { return fin::sqrt(q); }, x, y);
}
template <class T>
[[nodiscard]] inline C<T> sinh(T x, T y) noexcept {
    if (!std::isnan(x) && std::signbit(x)) return -sinh(-x, -y);
    if (negy(y)) return cj(sinh(x, -y));
    const T inf = std::numeric_limits<T>::infinity();
    if (fin(x)) {
        if (!fin(y)) return x == T() ? C<T>(x, y - y) : nn<T>();
        if (y == T()) return {detail::lm::sinh(x), y};
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
        if (y == T()) return {detail::lm::cosh(x), y};
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
    if (x == inf) return {T(1), fin(y) ? std::copysign(T(), detail::lm::sin(T(2) * y)) : T()};
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
// Spec grammar (the std-format-spec of floating-point types, split in two):
//   [[fill]align][sign]['#'][width][.precision]['L'][type]
// * fill, align and width apply to the whole "(a,b,c,d)" text; the default alignment is
//   right, as for numbers.  width and precision may be nested replacement fields ({} or {n}).
// * sign, '#', precision, 'L' and type (a A e E f F g G) apply to every component.
// * The '0' flag is rejected (format_error): zero padding of the whole text is meaningless.
//   std::format("{:*^30.2f}", q)  ->  "****(1.00,2.00,3.00,4.00)*****"
template <class T, class CharT>
struct std::formatter<quat::quaternion<T>, CharT> {
private:
    CharT fill_[4] = {CharT(' ')};
    int fill_len_ = 1;
    char align_ = 0;  // 0, '<', '>', '^'
    char sign_ = 0;   // 0, '+', '-', ' '
    bool alt_ = false, loc_ = false;
    int width_ = -1, width_arg_ = -1, prec_ = -1, prec_arg_ = -1;
    char type_ = 0;

    static constexpr bool is_align(CharT c) { return c == CharT('<') || c == CharT('>') || c == CharT('^'); }
    static constexpr int unit_len(const CharT* p, const CharT* e) {
        if constexpr (sizeof(CharT) == 1) {
            const unsigned char u = static_cast<unsigned char>(*p);
            const int n = u < 0x80 ? 1 : u >= 0xF0 ? 4 : u >= 0xE0 ? 3 : u >= 0xC0 ? 2 : 1;
            return n <= e - p ? n : 1;
        } else if constexpr (sizeof(CharT) == 2) {
            const auto u = static_cast<unsigned>(*p);
            return (u >= 0xD800 && u < 0xDC00 && e - p >= 2) ? 2 : 1;
        } else {
            return 1;
        }
    }
    template <class It>
    static constexpr It parse_int_or_arg(It it, It end, std::basic_format_parse_context<CharT>& ctx,
                                         int& value, int& arg) {
        if (*it == CharT('{')) {
            ++it;
            if (it == end) throw std::format_error("quaternion format: unterminated nested field");
            if (*it == CharT('}')) {
                arg = static_cast<int>(ctx.next_arg_id());
            } else {
                int n = 0;
                if (*it < CharT('0') || *it > CharT('9')) throw std::format_error("quaternion format: bad nested field");
                while (it != end && *it >= CharT('0') && *it <= CharT('9')) n = n * 10 + static_cast<int>(*it++ - CharT('0'));
                if (it == end || *it != CharT('}')) throw std::format_error("quaternion format: bad nested field");
                ctx.check_arg_id(static_cast<std::size_t>(n));
                arg = n;
            }
            return ++it;
        }
        int n = 0;
        while (it != end && *it >= CharT('0') && *it <= CharT('9')) {
            n = n * 10 + static_cast<int>(*it++ - CharT('0'));
            if (n > 1000000) throw std::format_error("quaternion format: number too large");
        }
        value = n;
        return it;
    }
    template <class Ctx>
    static int arg_value(Ctx& ctx, int id) {
        const auto a = ctx.arg(static_cast<std::size_t>(id));
        long long v = -1;
        std::visit_format_arg(
            [&v](auto x) {
                if constexpr (std::is_integral_v<decltype(x)> && !std::is_same_v<decltype(x), bool> &&
                              !std::is_same_v<decltype(x), CharT>)
                    v = (x < 0 || static_cast<unsigned long long>(x) > 1000000ULL) ? -1 : static_cast<long long>(x);
            },
            a);
        if (v < 0) throw std::format_error("quaternion format: width/precision argument must be a small non-negative integer");
        return static_cast<int>(v);
    }

public:
    constexpr auto parse(std::basic_format_parse_context<CharT>& ctx) {
        auto it = ctx.begin();
        const auto end = ctx.end();
        if (it == end || *it == CharT('}')) return it;
        // [[fill]align]
        const int ul = unit_len(&*it, &*it + (end - it));
        if (end - it > ul && is_align(it[ul])) {
            if (*it == CharT('{') || *it == CharT('}')) throw std::format_error("quaternion format: invalid fill");
            fill_len_ = ul;
            for (int n = 0; n < ul; ++n) fill_[n] = it[n];
            it += ul;
            align_ = static_cast<char>(*it++);
        } else if (is_align(*it)) {
            align_ = static_cast<char>(*it++);
        }
        if (it != end && (*it == CharT('+') || *it == CharT('-') || *it == CharT(' '))) sign_ = static_cast<char>(*it++);
        if (it != end && *it == CharT('#')) { alt_ = true; ++it; }
        if (it != end && *it == CharT('0'))
            throw std::format_error("quaternion format: the '0' flag is not supported; use a fill character");
        if (it != end && ((*it >= CharT('1') && *it <= CharT('9')) || *it == CharT('{')))
            it = parse_int_or_arg(it, end, ctx, width_, width_arg_);
        if (it != end && *it == CharT('.')) {
            ++it;
            if (it == end || !((*it >= CharT('0') && *it <= CharT('9')) || *it == CharT('{')))
                throw std::format_error("quaternion format: missing precision");
            it = parse_int_or_arg(it, end, ctx, prec_, prec_arg_);
        }
        if (it != end && *it == CharT('L')) { loc_ = true; ++it; }
        if (it != end && *it != CharT('}')) {
            switch (static_cast<char>(*it)) {
                case 'a': case 'A': case 'e': case 'E': case 'f': case 'F': case 'g': case 'G':
                    type_ = static_cast<char>(*it++);
                    break;
                default:
                    throw std::format_error("quaternion format: invalid type");
            }
        }
        if (it != end && *it != CharT('}')) throw std::format_error("quaternion format: invalid spec");
        return it;
    }

    template <class FormatContext>
    auto format(const quat::quaternion<T>& q, FormatContext& ctx) const {
        const int width = width_arg_ >= 0 ? arg_value(ctx, width_arg_) : width_;
        const int prec = prec_arg_ >= 0 ? arg_value(ctx, prec_arg_) : prec_;
        // Per-component spec "{:[sign][#][.prec][L][type]}".
        CharT spec[24];
        int n = 0;
        spec[n++] = CharT('{'); spec[n++] = CharT(':');
        if (sign_) spec[n++] = CharT(sign_);
        if (alt_) spec[n++] = CharT('#');
        if (prec >= 0) {
            spec[n++] = CharT('.');
            char d[12];
            int k = 0, p = prec;
            do { d[k++] = static_cast<char>('0' + p % 10); p /= 10; } while (p);
            while (k) spec[n++] = CharT(d[--k]);
        }
        if (loc_) spec[n++] = CharT('L');
        if (type_) spec[n++] = CharT(type_);
        spec[n++] = CharT('}');
        const std::basic_string_view<CharT> sv(spec, static_cast<std::size_t>(n));
        std::basic_string<CharT> body;
        body += CharT('(');
        const T c[4] = {q.w(), q.x(), q.y(), q.z()};
        for (int m = 0; m < 4; ++m) {
            if (m) body += CharT(',');
            const T v = c[m];
            if constexpr (std::is_same_v<CharT, char>) {
                if (loc_) std::vformat_to(std::back_inserter(body), ctx.locale(), sv, std::make_format_args(v));
                else std::vformat_to(std::back_inserter(body), sv, std::make_format_args(v));
            } else {
                if (loc_) std::vformat_to(std::back_inserter(body), ctx.locale(), sv, std::make_wformat_args(v));
                else std::vformat_to(std::back_inserter(body), sv, std::make_wformat_args(v));
            }
        }
        body += CharT(')');
        // Pad the whole text; its width counts code points (all components are ASCII
        // unless a locale supplies other separators).
        std::size_t len = 0;
        for (std::size_t m = 0; m < body.size(); m += static_cast<std::size_t>(unit_len(body.data() + m, body.data() + body.size()))) ++len;
        const std::size_t pad = width > 0 && static_cast<std::size_t>(width) > len ? static_cast<std::size_t>(width) - len : 0;
        const std::size_t before = align_ == '<' ? 0 : align_ == '^' ? pad / 2 : pad;
        auto out = ctx.out();
        const auto put_fill = [&](std::size_t k) {
            for (; k; --k)
                for (int f = 0; f < fill_len_; ++f) *out++ = fill_[f];
        };
        put_fill(before);
        for (const CharT ch : body) *out++ = ch;
        put_fill(pad - before);
        return out;
    }
};

#if defined(QUATERNION_MSVC_CONTRACT_OFF) && !defined(QUATERNION_DIAG_NORESTORE)
#pragma float_control(pop)
#endif

#endif  // QUATERNION_QUATERNION_HPP
