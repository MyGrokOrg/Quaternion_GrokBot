// Temporary macOS diagnostics (removed before merge).
#include <quaternion.hpp>
#include <cfloat>
#include <complex>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <limits>
#include <random>
using namespace quat;
template <class T> void sc(const char* name) {
    std::mt19937_64 g(7);
    long bad_s = 0, bad_l = 0, n = 0; int shown = 0;
    using L = std::numeric_limits<T>;
    for (int i = 0; i < 400000; ++i) {
        T x = std::ldexp(T(1) + T(g() % (1u << 30)) / T(1u << 30), int(g() % 80) - 40);
        if (g() & 1) x = -x;
        int e = int(g() % 120) + (L::min_exponent - L::digits - 60) + 40;  // push into subnormals
        if (i % 3 == 0) { x = std::ldexp(x, L::min_exponent - 20 + int(g() % 40)); e = int(g() % 200) - 100; }
        const T a = std::scalbn(x, e), b = detail::cx::scalbn(x, e);
        ++n;
        if (std::memcmp(&a, &b, sizeof(T) < 10 ? sizeof(T) : 10)) { ++bad_s; if (shown++ < 4) std::printf("  %s scalbn(%La, %d): libm %La, exact %La\n", name, (long double)x, e, (long double)a, (long double)b); }
        const T la = std::logb(x), lb = detail::cx::logb(x);
        if (la != lb) { ++bad_l; if (shown++ < 8) std::printf("  %s logb(%La): libm %La, cx %La\n", name, (long double)x, (long double)la, (long double)lb); }
    }
    std::printf("%s: scalbn differs in %ld of %ld, logb in %ld\n", name, bad_s, n, bad_l);
}
template <class T> void inv(const char* name) {
    const T pts[][2] = {{T(0.25), T(0.5)}, {T(-0.75), T(1.5)}, {T(2.5), T(-0.25)}, {T(-3), T(-2)}, {T(0.5), T(10)}, {T(1e-3), T(2e-3)}};
    const quaternion<T> u = quaternion<T>(T(0), T(2) / 3, T(-1) / 3, T(2) / 3);
    for (auto& p : pts) {
        const std::complex<T> z(p[0], p[1]);
#define D(fn) { const std::complex<T> w = std::fn(z); const quaternion<T> q = quat::fn(quaternion<T>(p[0]) + p[1] * u); \
        const quaternion<T> want = quaternion<T>(w.real()) + w.imag() * u; const double err = double(abs(q - want) / abs(want)) / double(std::numeric_limits<T>::epsilon()); \
        if (err > 8) std::printf("  %s %s(%g%+gi): %.3g eps  quat %.17Lg %+.17Lg  std::complex %.17Lg %+.17Lg\n", name, #fn, double(p[0]), double(p[1]), err, (long double)q.real(), (long double)abs(q.unreal()), (long double)w.real(), (long double)w.imag()); }
        D(asin) D(acos) D(atan) D(asinh) D(acosh) D(atanh)
#undef D
    }
}
int main() {
    sc<double>("double"); sc<long double>("long double"); sc<float>("float");
    inv<double>("double"); inv<long double>("long double");
    std::printf("sizeof(long double) %zu, digits %d, FLT_EVAL %d\n", sizeof(long double), std::numeric_limits<long double>::digits, (int)FLT_EVAL_METHOD);
}
