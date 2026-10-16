// The library's own x87 long double kernels (QUATERNION_X87_LIBM, quaternion/detail/x87_libm.hpp).
// Built everywhere; the checks run only where the macro is on (CI forces it on Linux x86_64
// with -DQUATERNION_X87_LIBM=1, and it is the default on macOS x86_64).
//  * every function against correctly rounded mpmath references (tests/x87_refs.inc): the
//    result must be the reference or a neighbour (error below 1 ulp; the kernels measure about 0.5), and the share that is
//    exactly the reference is printed;
//  * special values (+-0, +-inf, NaN) behave as the C library's;
//  * where the C library is glibc, a seeded random comparison with it (glibc's x87 long double
//    functions are accurate to a few ulp; the bounds are 3 ulp, 8 for tgamma and lgamma);
//  * the overloads in quat::detail::lm are the ones the library calls.
#include "test_framework.hpp"

#if QUATERNION_X87_LIBM
#include <cstdint>
#include <cstring>
#include <random>
#include <string_view>

namespace {
using ld = long double;
namespace lm = quat::detail::lm;

struct ref {
    const char* f;
    ld x, y, r;
};
constexpr ld inf = std::numeric_limits<ld>::infinity();
const ref refs[] = {
#include "x87_refs.inc"
};

ld call(std::string_view f, ld x, [[maybe_unused]] ld y) {
    if (f == "sin") return lm::sin(x);
    if (f == "cos") return lm::cos(x);
    if (f == "tan") return lm::tan(x);
    if (f == "sinh") return lm::sinh(x);
    if (f == "expm1") return lm::expm1(x);
    if (f == "pow") return lm::pow(x, y);
    if (f == "tgamma") return lm::tgamma(x);
    if (f == "lgamma") return lm::lgamma(x);
    return std::numeric_limits<ld>::quiet_NaN();
}
ld call_std(std::string_view f, ld x, [[maybe_unused]] ld y) {
    if (f == "sin") return std::sin(x);
    if (f == "cos") return std::cos(x);
    if (f == "tan") return std::tan(x);
    if (f == "sinh") return std::sinh(x);
    if (f == "expm1") return std::expm1(x);
    if (f == "pow") return std::pow(x, y);
    if (f == "tgamma") return std::tgamma(x);
    if (f == "lgamma") return std::lgamma(x);
    return std::numeric_limits<ld>::quiet_NaN();
}
// distance in units of the last place of b (b finite, nonzero)
ld ulps(ld a, ld b) {
    int e;
    std::frexp(b, &e);
    const ld u = std::ldexp(1.0L, (e - 64 > -16445 ? e - 64 : -16445));
    return std::fabs(a - b) / u;
}
bool same_bits(ld a, ld b) {
    unsigned char p[10], q[10];
    std::memcpy(p, &a, 10);
    std::memcpy(q, &b, 10);
    return std::memcmp(p, q, 10) == 0;
}
}  // namespace

TEST(mpmath_references) {
    struct st { std::string_view f; int n = 0, exact = 0; ld worst = 0; };
    st s[8];
    int ns = 0;
    for (const ref& t : refs) {
        const ld got = call(t.f, t.x, t.y);
        const ld d = t.r == 0 || std::isinf(t.r) ? (same_bits(got, t.r) ? 0.0L : inf) : ulps(got, t.r);
        if (!(d <= 1.0L))
            std::fprintf(stderr, "  %s(%La, %La) = %La, reference %La (%Lg ulp)\n", t.f, t.x, t.y, got, t.r, d);
        CHECK(d <= 1.0L);
        int i = 0;
        while (i < ns && s[i].f != t.f) ++i;
        if (i == ns) s[ns++].f = t.f;
        ++s[i].n;
        s[i].exact += same_bits(got, t.r);
        s[i].worst = std::max(s[i].worst, d);
    }
    for (int i = 0; i < ns; ++i)
        std::printf("  %-7.*s %4d references, %4d correctly rounded, the others 1 ulp off\n",
                    static_cast<int>(s[i].f.size()), s[i].f.data(), s[i].n, s[i].exact);
}

TEST(special_values) {
    const ld nan = std::numeric_limits<ld>::quiet_NaN(), big = std::numeric_limits<ld>::max();
    const ld tiny = std::numeric_limits<ld>::denorm_min(), mn = std::numeric_limits<ld>::min();
    const ld xs[] = {0.0L, -0.0L, inf, -inf, nan, tiny, -tiny, mn, -mn, big, -big, 1.0L, -1.0L, 2.0L, -2.0L, 0.5L,
                     -0.5L, 11356.0L, 11357.0L, -11400.0L, 1754.0L, 1755.5L, -1754.5L, -1760.5L, 0x1p-70L, -0x1p-70L};
    auto agree = [](ld a, ld b) {
        return (std::isnan(a) && std::isnan(b)) || same_bits(a, b) ||
               (std::isfinite(b) && b != 0 && std::isfinite(a) && ulps(a, b) <= 2.0L);
    };
    for (const char* f : {"sin", "cos", "tan", "sinh", "expm1", "tgamma", "lgamma"})
        for (ld x : xs) {
            const ld a = call(f, x, 0), b = call_std(f, x, 0);
            if (!agree(a, b)) std::fprintf(stderr, "  %s(%La) = %La, C library %La\n", f, x, a, b);
            CHECK(agree(a, b));
        }
    for (ld x : xs)
        for (ld y : xs) {
            const ld a = lm::pow(x, y), b = std::pow(x, y);
            if (!agree(a, b)) std::fprintf(stderr, "  pow(%La, %La) = %La, C library %La\n", x, y, a, b);
            CHECK(agree(a, b));
        }
    CHECK(std::signbit(lm::sin(-0.0L)) && std::signbit(lm::tan(-0.0L)) && lm::cos(-0.0L) == 1.0L);
    CHECK(std::signbit(lm::sinh(-0.0L)) && std::signbit(lm::expm1(-0.0L)));
    CHECK(lm::pow(-2.0L, 3) == -8.0L && lm::pow(2, 0.5L) == lm::pow(2.0L, 0.5L));
}

TEST(overloads_are_used) {
    // the non-template long double overloads win over the forwarders, and float/double still
    // go to std::
    static_assert(std::is_same_v<decltype(lm::sin(1.0L)), long double>);
    const ld x = 0x1.921fb54442d18469p+1L * 1e6L;  // pi * 1e6, rounded
    CHECK(same_bits(lm::sin(x), quat::detail::lm::x87::trig_eval(x, quat::detail::lm::x87::trig::sin)));
    CHECK(lm::sin(1.0) == std::sin(1.0) && lm::cos(1.0f) == std::cos(1.0f));
    // ... and the library's functions go through them: sin(a + 0 i) for a quaternion
    const quat::quaternion<ld> q(x, 0, 0, 0);
    CHECK(same_bits(quat::sin(q).w(), lm::sin(x)));
}

#if defined(__GLIBC__)
TEST(random_vs_glibc) {
    // seed 20261016; 20,000 inputs per function (pow: 20,000 pairs), drawn below.  The kernels
    // are within about 0.55 ulp of mpmath, so the differences are glibc's own errors: up to
    // about 1.4 ulp for the trigonometric functions, sinh, expm1 and pow, and several ulp for
    // tgamma and lgamma; the bounds are 3 and 8 ulp.  lgamma for x < 0 is compared only where
    // |lgamma| > 1/2 (glibc loses relative accuracy next to its zeros there).
    std::mt19937_64 g(20261016);
    auto u = [&g](ld a, ld b) { return a + (b - a) * static_cast<ld>(g() >> 11) * 0x1p-53L; };
    auto mag = [&g](int lo, int hi) {
        const int e = lo + static_cast<int>(g() % static_cast<std::uint64_t>(hi - lo));
        const ld x = std::ldexp(static_cast<ld>(g() | (std::uint64_t(1) << 63)), e - 64);
        return (g() & 1) ? -x : x;
    };
    for (const char* f : {"sin", "cos", "tan", "sinh", "expm1", "pow", "tgamma", "lgamma"}) {
        const std::string_view fs = f;
        ld worst = 0;
        int differ = 0, n = 0;
        for (int i = 0; i < 20000; ++i) {
            ld x = 0, y = 0;
            const bool lo = i < 10000;
            if (fs == "sin" || fs == "cos" || fs == "tan") x = lo ? mag(-40, 16) : mag(16, 16384);
            else if (fs == "sinh") x = lo ? mag(-40, 3) : u(-11350, 11350);
            else if (fs == "expm1") x = lo ? mag(-40, 3) : u(-60, 11350);
            else if (fs == "pow") {
                x = lo ? u(0.01L, 100) : std::fabs(mag(-60, 60));
                y = lo ? u(-60, 60) : u(-150, 150);
            } else if (fs == "tgamma") x = lo ? u(0.001L, 1750) : u(-170, 0);
            else x = lo ? u(0.001L, 1000) : u(-50, 0);
            const ld a = call(f, x, y), b = call_std(f, x, y);
            if (!std::isfinite(b) || b == 0 || (fs == "lgamma" && x < 0 && std::fabs(b) < 0.5L)) continue;
            ++n;
            const ld d = ulps(a, b);
            differ += !same_bits(a, b);
            worst = std::max(worst, d);
        }
        const ld bound = (fs == "tgamma" || fs == "lgamma") ? 8.0L : 3.0L;
        std::printf("  %-6s vs glibc: %5d of %5d differ, at most %Lg ulp (bound %Lg)\n", f, differ, n, worst, bound);
        CHECK(worst <= bound);
    }
}
#endif

#else
TEST(skipped) { std::printf("  QUATERNION_X87_LIBM is 0 in this build: nothing to test\n"); }
#endif

TEST_MAIN
