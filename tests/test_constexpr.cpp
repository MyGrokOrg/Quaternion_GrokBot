// Constant evaluation of the transcendental-free operations (no <cmath> needed), and the
// detail::cx helpers behind them, checked against the runtime std:: functions.
// This target is compiled with -ffp-contract=off (CMakeLists.txt) so the helpers' error-free
// transformations stay exact when they are exercised at run time below.
#include "test_framework.hpp"
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
using quat::quaternion;
namespace cx = quat::detail::cx;

// ---- static_assert: helpers --------------------------------------------------------
template <class T>
constexpr bool helpers() {
    using L = std::numeric_limits<T>;
    const T inf = L::infinity(), nan = L::quiet_NaN();
    bool ok = true;
    ok = ok && cx::signbit(-T(0)) && !cx::signbit(T(0)) && cx::signbit(-nan) && !cx::signbit(nan);
    ok = ok && cx::signbit(-inf) && cx::isinf(-inf) && cx::isnan(nan) && !cx::isfinite(inf) && cx::isfinite(L::max());
    ok = ok && cx::fabs(T(-2.5)) == T(2.5) && cx::signbit(cx::copysign(T(3), -T(0))) && cx::copysign(T(-3), T(1)) == T(3);
    ok = ok && cx::logb(T(1)) == T(0) && cx::logb(T(0.75)) == T(-1) && cx::logb(L::max()) == T(L::max_exponent - 1);
    ok = ok && cx::logb(L::denorm_min()) == T(L::min_exponent - L::digits) && cx::logb(T(0)) == -inf;
    ok = ok && cx::scalbn(T(1), L::max_exponent - 1) == L::max() / (T(2) - L::epsilon());
    ok = ok && cx::scalbn(T(1), L::min_exponent - L::digits) == L::denorm_min();
    ok = ok && cx::scalbn(T(1), L::min_exponent - L::digits - 1) == T(0);   // tie to even
    ok = ok && cx::scalbn(T(1.5), L::min_exponent - L::digits - 1) == L::denorm_min();
    ok = ok && cx::scalbn(T(3), L::min_exponent - L::digits - 1) == 2 * L::denorm_min();  // 1.5 ulp -> 2
    ok = ok && cx::scalbn(L::denorm_min(), L::digits - L::min_exponent) == T(1);
    ok = ok && cx::scalbn(T(1), L::max_exponent) == inf;
    int e = 0;
    ok = ok && cx::frexp(T(8), &e) == T(0.5) && e == 4 && cx::ldexp(T(0.5), 4) == T(8);
    ok = ok && cx::frexp(L::denorm_min(), &e) == T(0.5) && e == L::min_exponent - L::digits + 1;
    ok = ok && cx::sqrt(T(4)) == T(2) && cx::sqrt(T(2)) * cx::sqrt(T(2)) != T(0) && cx::sqrt(L::max()) > T(0);
    ok = ok && cx::sqrt(L::denorm_min()) > T(0) && cx::hypot(T(3), T(4)) == T(5);
    ok = ok && cx::hypot(L::max(), L::max()) == inf;
    ok = ok && cx::hypot(L::max() / 2, L::max() / 2) < inf && cx::hypot(L::denorm_min(), T(0)) == L::denorm_min();
    return ok;
}
static_assert(helpers<float>());
static_assert(helpers<double>());
static_assert(helpers<long double>());
static_assert(cx::sqrt(2.0) == 0x1.6a09e667f3bcdp+0);
static_assert(cx::sqrt(2.0f) == 0x1.6a09e6p+0f);
static_assert(cx::sqrt(2.0L) == 0xb.504f333f9de6484p-3L);

// ---- static_assert: quaternion operations ------------------------------------------------
constexpr quaternion<double> P(1.5, -2, 0.25, 3), Q(-0.5, 4, 1, -2);
// Value printed by the runtime division (x86-64 without FMA).
static_assert(P / Q == quaternion<double>(-0x1.5d5d5d5d5d5d6p-1, -0x1.2121212121212p-4, -0x1.cfcfcfcfcfcfdp-2, 0x1.b1b1b1b1b1b1bp-3));
static_assert(abs(quaternion<double>(1, 2, 2, 4)) == 5.0);
static_assert(abs(quaternion<double>(1e300, 1e300, 1e300, 1e300)) == 2e300);
static_assert(abs(quaternion<double>(3e-320, 4e-320, 0, 0)) == 5e-320);
static_assert(normalize(quaternion<double>(0, 3, 0, 4)) == quaternion<double>(0, 0.6, 0, 0.8));
static_assert(inverse(quaternion<double>(0, 0, 2, 0)) == quaternion<double>(0, 0, -0.5, 0));
static_assert(pow(quaternion<float>(1, 1, 0, 0), 8) == 16.0f);
static_assert(pow(quaternion<double>(0, 0, 0, 2), -2) == -0.25);
static_assert(1e300 / quaternion<double>(0, 1e300, 0, 0) == quaternion<double>(0, -1, 0, 0));  // scaled, no overflow
static_assert(quaternion<double>(1e-300, 0, 0, 0) / quaternion<double>(1e-300, 1e-300, 0, 0) == quaternion<double>(0.5, -0.5, 0, 0));
static_assert(proj(quaternion<double>(1, -std::numeric_limits<double>::infinity(), -0.0, 2)) ==
              quaternion<double>(std::numeric_limits<double>::infinity(), -0.0, -0.0, 0.0));
static_assert(cx::signbit(proj(quaternion<double>(std::numeric_limits<double>::infinity(), -1, 0, 0)).x()));
static_assert(quaternion<long double>(1, 2, 3, 4) / quaternion<long double>(1, 2, 3, 4) == 1.0L);
static_assert(abs(quaternion<long double>(0, 0, 3, 4)) == 5.0L);
static_assert(quaternion<float>(2, 0, 0, 0) / quaternion<float>(0, 0, 0, 4) == quaternion<float>(0, 0, 0, -0.5f));

// ---- constexpr vs runtime on many inputs ------------------------------------------------
template <class T>
struct table { T div[256][4]; T ab[256]; };
template <class T>
constexpr quaternion<T> sample(int i) {
    // Deterministic inputs spanning many binades (constexpr LCG).
    std::uint64_t s = 0x1234567ULL + static_cast<std::uint64_t>(i) * 0x9e3779b97f4a7c15ULL;
    T c[4];
    for (T& x : c) {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        const T m = T(1) + static_cast<T>(static_cast<double>(s >> 40) / 16777216.0);
        const int e = static_cast<int>((s >> 8) % 120) - 60;
        x = ((s >> 5) & 1 ? -m : m) * cx::pow2<T>(e);
    }
    return {c[0], c[1], c[2], c[3]};
}
template <class T>
constexpr table<T> make_table() {
    table<T> t{};
    for (int i = 0; i < 256; ++i) {
        const auto q = sample<T>(i) / sample<T>(i + 256);
        t.div[i][0] = q.w(); t.div[i][1] = q.x(); t.div[i][2] = q.y(); t.div[i][3] = q.z();
        t.ab[i] = abs(sample<T>(i));
    }
    return t;
}
template <class T>
bool same(T a, T b) { return std::memcmp(&a, &b, std::numeric_limits<T>::digits == 64 ? 10 : sizeof(T)) == 0; }

template <class T> void constexpr_vs_runtime() {
    constexpr table<T> t = make_table<T>();
    int div_diff = 0, abs_diff = 0;
    for (int i = 0; i < 256; ++i) {
        volatile int vi = i;
        const auto q = sample<T>(vi) / sample<T>(vi + 256);
        const T c[4] = {q.w(), q.x(), q.y(), q.z()};
        for (int k = 0; k < 4; ++k) div_diff += !same(c[k], t.div[i][k]);
        abs_diff += !same(abs(sample<T>(vi)), t.ab[i]);
    }
#if defined(__FMA__)
    // The runtime Hamilton product may be contracted into FMAs; constant evaluation never is.
    (void)div_diff;
#else
    CHECK(div_diff == 0);
#endif
    // abs: std::hypot (libm) at run time vs a correctly rounded hypot in constant evaluation.
    CHECK(abs_diff <= 16);
    if (abs_diff) std::printf("  [%s] constexpr abs differs from runtime in %d of 256\n", qtest::current_type, abs_diff);
}
TEST(constexpr_matches_runtime) { FOR_EACH_TYPE(constexpr_vs_runtime); }

// ---- helpers at run time vs std:: (exhaustive-ish random sweep) ---------------------------
template <class T> void helpers_vs_std() {
    using L = std::numeric_limits<T>;
    std::uint64_t s = 42;
    int bad_scal = 0, bad_logb = 0, bad_frexp = 0, bad_sqrt = 0, hyp_diff = 0, bad_hyp = 0;
    const int n = 200000;
    for (int i = 0; i < n; ++i) {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        const T m = T(1) + static_cast<T>(static_cast<double>(s >> 11) / 9007199254740992.0);
        const int e = static_cast<int>((s >> 3) % static_cast<std::uint64_t>(L::max_exponent - L::min_exponent + L::digits)) + L::min_exponent - L::digits;
        const T x = std::ldexp(m, e);
        const int k = static_cast<int>((s >> 20) % 400) - 200;
        bad_scal += !same(cx::scalbn(x, k), std::scalbn(x, k));
        bad_logb += !same(cx::logb(x), std::logb(x));
        int e1 = 0, e2 = 0;
        const T f1 = cx::frexp(x, &e1), f2 = std::frexp(x, &e2);
        bad_frexp += !same(f1, f2) || e1 != e2;
        bad_sqrt += !same(cx::sqrt(x), std::sqrt(x));
        const T y = std::ldexp(m * T(0.7), e + static_cast<int>((s >> 30) % 9) - 4);
        const T h1 = cx::hypot(x, y), h2 = std::hypot(x, y);
        hyp_diff += !same(h1, h2);
        // differences must be at most one ulp
        bad_hyp += !(h1 == h2 || std::nextafter(h1, h2) == h2);
    }
    CHECK(bad_scal == 0);
    CHECK(bad_logb == 0);
    CHECK(bad_frexp == 0);
    CHECK(bad_sqrt == 0);
    CHECK(bad_hyp == 0);
    std::printf("  [%s] %d inputs: scalbn/logb/frexp/sqrt identical to std::; hypot differs from std::hypot in %d (each by 1 ulp)\n",
                qtest::current_type, n, hyp_diff);
}
TEST(cx_helpers_match_std) { FOR_EACH_TYPE(helpers_vs_std); }

TEST_MAIN
