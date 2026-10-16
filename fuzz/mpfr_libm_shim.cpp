// Informational (.github/workflows/libm-check.yml; never part of the library): correctly
// rounded replacements (MPFR, 64-bit precision) for exactly the long double C library functions that are inaccurate on macOS x86_64 (sinl, cosl,
// tanl, sinhl, expm1l, powl, tgammal, lgammal).  Linked into the executable, they take
// precedence over libm, so the library and libc++'s std::complex<long double> use them.
// Shows what the library reaches on macOS x86_64 once only these functions are accurate.
// Needs MPFR (LGPL), so it is linked into test executables only, never into the headers.
#include <mpfr.h>
namespace {
struct M {
    mpfr_t r, a, b;
    M() { mpfr_inits2(64, r, a, b, (mpfr_ptr)0); }
    ~M() { mpfr_clears(r, a, b, (mpfr_ptr)0); }
};
template <class F> long double un(long double x, F f) {
    M m; mpfr_set_ld(m.a, x, MPFR_RNDN);
    int t = f(m.r, m.a, MPFR_RNDN); mpfr_subnormalize(m.r, t, MPFR_RNDN);
    return mpfr_get_ld(m.r, MPFR_RNDN);
}
}  // namespace
extern "C" {
long double sinl(long double x) { return un(x, mpfr_sin); }
long double cosl(long double x) { return un(x, mpfr_cos); }
long double tanl(long double x) { return un(x, mpfr_tan); }
long double sinhl(long double x) { return un(x, mpfr_sinh); }
long double expm1l(long double x) { return un(x, mpfr_expm1); }
long double tgammal(long double x) { return un(x, mpfr_gamma); }
long double lgammal(long double x) {
    return un(x, [](mpfr_ptr r, mpfr_srcptr a, mpfr_rnd_t d) { int s; return mpfr_lgamma(r, &s, a, d); });
}
long double powl(long double x, long double y) {
    M m; mpfr_set_ld(m.a, x, MPFR_RNDN); mpfr_set_ld(m.b, y, MPFR_RNDN);
    int t = mpfr_pow(m.r, m.a, m.b, MPFR_RNDN); mpfr_subnormalize(m.r, t, MPFR_RNDN);
    return mpfr_get_ld(m.r, MPFR_RNDN);
}
}
