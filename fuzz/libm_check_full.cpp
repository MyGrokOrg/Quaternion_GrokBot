// Every long double C library function the headers call whose result is rounded (and fma), on
// several input ranges, printed in hex for fuzz/libm_check_full.py (mpmath).  Informational, not
// part of CI's test suites; .github/workflows/libm-check.yml runs it on demand (macOS x86_64 and
// Linux).  Inputs are built without libm calls, so every platform gets the same inputs.
//   c++ -std=c++23 -O2 fuzz/libm_check_full.cpp -o libm_check_full
//   ./libm_check_full 12 300 > out.txt && python3 fuzz/libm_check_full.py out.txt   (seed, samples)
//   With -DLIBM_CHECK_LM -Iinclude, the functions the library replaces under
//   QUATERNION_X87_LIBM (sin, cos, tan, sinh, expm1, pow, tgamma, lgamma; asin and acos for
//   |x| > 1/2) are called through
//   quat::detail::lm instead, i.e. what the headers use (the library's own kernels where the
//   macro is on, the C library's otherwise).
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#ifdef LIBM_CHECK_LM
#include <quaternion/quaternion.hpp>
namespace F = quat::detail::lm;
#else
namespace F = std;
#endif
int main(int argc, char** argv) {
    const unsigned long long seed = argc > 1 ? std::strtoull(argv[1], nullptr, 10) : 12;
    const int n = argc > 2 ? std::atoi(argv[2]) : 300;
    std::mt19937_64 g(seed);
    // full 64-bit significand in [a, b)
    auto U = [&](long double a, long double b) {
        long double r = (long double)(g() >> 1) / 9223372036854775808.0L;
        return a + (b - a) * r;
    };
    auto LU = [&](int lo2, int hi2) {  // sign * [1, 2) * 2^k, k uniform in [lo2, hi2), exact (no libm)
        long double s = (g() & 1) ? -1.0L : 1.0L;
        int k = lo2 + (int)(g() % (unsigned long long)(hi2 - lo2));
        return s * std::ldexp(U(1, 2), k);
    };
    auto P1 = [](const char* f, long double x, long double r) { std::printf("%s %La %La\n", f, x, r); };
    auto P2 = [](const char* f, long double x, long double y, long double r) {
        std::printf("%s %La %La %La\n", f, x, y, r);
    };
    std::printf("# seed %llu n %d LDBL_MANT_DIG %d\n", seed, n, (int)__LDBL_MANT_DIG__);
#ifdef LIBM_CHECK_LM
    std::printf("# through quat::detail::lm, QUATERNION_X87_LIBM %d\n", (int)QUATERNION_X87_LIBM);
#endif
    for (int i = 0; i < n; ++i) {
        long double x;
        x = U(-4, 4);            P1("sin@small", x, F::sin(x)); P1("cos@small", x, F::cos(x)); P1("tan@small", x, F::tan(x));
        x = U(-700, 700);        P1("sin@mid", x, F::sin(x));   P1("cos@mid", x, F::cos(x));   P1("tan@mid", x, F::tan(x));
        x = LU(10, 62);          P1("sin@large", x, F::sin(x)); P1("cos@large", x, F::cos(x)); P1("tan@large", x, F::tan(x));
        x = U(-40, 40);          P1("exp@mid", x, std::exp(x));
        x = U(-11000, 11000);    P1("exp@wide", x, std::exp(x));
        x = U(-40, 40);          P1("exp2@mid", x, std::exp2(x));
        x = U(-0.05L, 0.05L);    P1("expm1@small", x, F::expm1(x));
        x = U(-40, 40);          P1("expm1@mid", x, F::expm1(x));
        x = U(0.5L, 2);          P1("log@near1", x, std::log(x));
        x = std::fabs(LU(-16000, 16000)); P1("log@wide", x, std::log(x));
        x = U(-0.05L, 0.05L);    P1("log1p@small", x, std::log1p(x));
        x = U(-0.9L, 40);        P1("log1p@mid", x, std::log1p(x));
        x = U(-1, 1);            P1("sinh@small", x, F::sinh(x)); P1("cosh@small", x, std::cosh(x)); P1("tanh@small", x, std::tanh(x));
        x = U(-40, 40);          P1("sinh@mid", x, F::sinh(x)); P1("cosh@mid", x, std::cosh(x)); P1("tanh@mid", x, std::tanh(x));
        x = U(-11000, 11000);    P1("sinh@wide", x, F::sinh(x)); P1("cosh@wide", x, std::cosh(x));
        x = U(-1, 1);            P1("asin", x, F::asin(x)); P1("acos", x, F::acos(x));
        x = U(-5, 5);            P1("atan@mid", x, std::atan(x));
        x = LU(-60, 60);         P1("atan@wide", x, std::atan(x));
        { long double a = LU(-30, 30), b = LU(-30, 30); P2("atan2", a, b, std::atan2(a, b)); }
        { long double a = U(-40, 40), b = U(-5, 5);     P2("atan2@mid", a, b, std::atan2(a, b)); }
        { long double b = U(0.5L, 3), e = U(-40, 40);   P2("pow@mid", b, e, F::pow(b, e)); }
        { long double b = (long double)(2 + g() % 999), e = U(-60, 0); P2("pow@intbase", b, e, F::pow(b, e)); }
        { long double b = std::fabs(LU(-20, 20)), e = U(-500, 500); P2("pow@wide", b, e, F::pow(b, e)); }
        x = U(0.01L, 60);        P1("lgamma@mid", x, F::lgamma(x));
        x = U(60, 10000);        P1("lgamma@large", x, F::lgamma(x));
        x = U(-50, 0);           P1("lgamma@neg", x, F::lgamma(x));
        x = U(2.5L, 15);         P1("tgamma@mid", x, F::tgamma(x));
        x = U(0.01L, 170);       P1("tgamma@wide", x, F::tgamma(x));
        x = U(-50, 0);           P1("tgamma@neg", x, F::tgamma(x));
        x = U(-40, 40);          P1("cbrt", x, std::cbrt(x));
        x = U(0.01L, 60);        P1("sqrt", x, std::sqrt(x));
        { long double a = LU(-100, 100), b = a * LU(-30, 30); P2("hypot", a, b, std::hypot(a, b)); }
        {   // fma: random, and with cancellation c ~ -a*b (as two_prod uses it)
            long double a = LU(-20, 20), b = LU(-20, 20), c = LU(-20, 20);
            std::printf("fma@rand %La %La %La %La\n", a, b, c, std::fma(a, b, c));
            long double p = a * b;
            std::printf("fma@twoprod %La %La %La %La\n", a, b, -p, std::fma(a, b, -p));
        }
    }
    // asin and acos next to +-1: x = +-(1 - d), d = [1, 2) 2^-k, k uniform in [2, 41) (the
    // nightly fuzz caught Apple's asinl at 1 - |x| from 6e-12 to 9e-7).  Drawn after the loop
    // above, so the inputs of every other range are those of earlier versions.
    for (int i = 0; i < n; ++i) {
        const long double d = std::ldexp(U(1, 2), -2 - (int)(g() % 39ULL));
        const long double x = (g() & 1) ? -(1 - d) : 1 - d;
        P1("asin@pm1", x, F::asin(x)); P1("acos@pm1", x, F::acos(x));
    }
}
