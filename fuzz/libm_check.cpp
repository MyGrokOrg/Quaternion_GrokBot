// Long double C library functions on random inputs, printed in hex for fuzz/libm_check.py (mpmath).
// Informational, not part of CI: it measured Apple's x87 long double libm for the README (macOS x86_64).
#include <cmath>
#include <cstdio>
#include <random>
int main() {
    std::mt19937_64 g(11);
    auto U = [&](long double a, long double b) { return a + (b - a) * (long double)(g() >> 11) / 9007199254740992.0L; };
    for (int i = 0; i < 300; ++i) {
        long double x = U(-40, 40), y = U(0.01L, 60), s = U(-700, 700), t = U(0.5L, 3), l = U(-5, 5);
        std::printf("exp %La %La\n", x, std::exp(x));
        std::printf("expm1 %La %La\n", l * 0.01L, std::expm1(l * 0.01L));
        std::printf("log %La %La\n", y, std::log(y));
        std::printf("log1p %La %La\n", l * 0.01L, std::log1p(l * 0.01L));
        std::printf("sin %La %La\n", s, std::sin(s));
        std::printf("cos %La %La\n", s, std::cos(s));
        std::printf("sinh %La %La\n", x, std::sinh(x));
        std::printf("cosh %La %La\n", x, std::cosh(x));
        std::printf("atan2 %La %La %La\n", x, l, std::atan2(x, l));
        std::printf("pow %La %La %La\n", t, x, std::pow(t, x));
        std::printf("lgamma %La %La\n", y, std::lgamma(y));
        std::printf("tgamma %La %La\n", t * 5, std::tgamma(t * 5));
        std::printf("sqrt %La %La\n", y, std::sqrt(y));
        std::printf("cbrt %La %La\n", x, std::cbrt(x));
        std::printf("atan %La %La\n", l, std::atan(l));
    }
}
