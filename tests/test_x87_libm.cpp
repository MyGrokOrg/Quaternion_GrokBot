// The library's own x87 long double kernels (QUATERNION_X87_LIBM, quaternion/detail/x87_libm.hpp).
// Built everywhere; the checks run only where the macro is on (CI forces it on Linux x86_64
// with -DQUATERNION_X87_LIBM=1, and it is the default on macOS x86_64).
//  * every function against correctly rounded mpmath references (tests/x87_refs.inc): the
//    result must be the reference or a neighbour (error below 1 ulp; the kernels measure about 0.5), and the share that is
//    exactly the reference is printed;
//  * special values (+-0, +-inf, NaN) behave as the C library's;
//  * where the C library is glibc, a seeded random comparison with it (glibc's x87 long double
//    functions are accurate to about 1 ulp, so the two must stay within 2 ulp of each other);
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
    return std::numeric_limits<ld>::quiet_NaN();
}
ld call_std(std::string_view f, ld x, [[maybe_unused]] ld y) {
    if (f == "sin") return std::sin(x);
    if (f == "cos") return std::cos(x);
    if (f == "tan") return std::tan(x);
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
    const ld nan = std::numeric_limits<ld>::quiet_NaN();
    for (const char* f : {"sin", "cos", "tan"}) {
        for (ld x : {0.0L, -0.0L, inf, -inf, nan, std::numeric_limits<ld>::denorm_min(), -std::numeric_limits<ld>::min(),
                     std::numeric_limits<ld>::max(), -std::numeric_limits<ld>::max()}) {
            const ld a = call(f, x, 0), b = call_std(f, x, 0);
            CHECK((std::isnan(a) && std::isnan(b)) || same_bits(a, b) || ulps(a, b) <= 1.0L);
        }
    }
    CHECK(std::signbit(lm::sin(-0.0L)) && std::signbit(lm::tan(-0.0L)) && lm::cos(-0.0L) == 1.0L);
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
    // seed 20261016; 20,000 inputs per function: 10,000 with |x| in [2^-40, 2^16) and 10,000 in
    // [2^16, 2^16384), random significands and signs
    std::mt19937_64 g(20261016);
    for (const char* f : {"sin", "cos", "tan"}) {
        ld worst = 0;
        int differ = 0;
        for (int i = 0; i < 20000; ++i) {
            const int e = i < 10000 ? static_cast<int>(g() % 56) - 40 : 16 + static_cast<int>(g() % 16368);
            ld x = std::ldexp(static_cast<ld>(g() | (std::uint64_t(1) << 63)), e - 64);
            if (g() & 1) x = -x;
            const ld a = call(f, x, 0), b = call_std(f, x, 0);
            const ld d = ulps(a, b);
            differ += !same_bits(a, b);
            worst = std::max(worst, d);
        }
        std::printf("  %-4s vs glibc: %d of 20000 differ, at most %Lg ulp\n", f, differ, worst);
        CHECK(worst <= 2.0L);
    }
}
#endif

#else
TEST(skipped) { std::printf("  QUATERNION_X87_LIBM is 0 in this build: nothing to test\n"); }
#endif

TEST_MAIN
