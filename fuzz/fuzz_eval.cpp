// Evaluator for the nightly mpmath fuzz job (fuzz/fuzz.py).  Reads lines
//   fn p1 p2 a b axis
// (hex floats; axis 1..3 = i, j, k, negative = reversed; 4 = the oblique unit vector
// u = (2, 3, 6)/7, i.e. vector part (2c, 3c, 6c) with c = b/7, exact when b = 7c with c
// short enough) and evaluates fn on the quaternion a + b u in the type given by argv[1]
// ("double" or "long double").  Prints the real part and the component along u in hex, or
// "nan nan" for an unknown function.
#include <quaternion.hpp>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>

template <class T>
static quat::quaternion<T> eval(const std::string& f, T p1, T p2, const quat::quaternion<T>& q) {
    namespace Q = quat;
    const int n1 = static_cast<int>(p1);
    const unsigned u1 = static_cast<unsigned>(p1);
    // elementary functions (quaternion.hpp)
    if (f == "exp") return Q::exp(q);
    if (f == "log") return Q::log(q);
    if (f == "log10") return Q::log10(q);
    if (f == "sqrt") return Q::sqrt(q);
    if (f == "pow") return Q::pow(q, p1);
    if (f == "sin") return Q::sin(q);
    if (f == "cos") return Q::cos(q);
    if (f == "tan") return Q::tan(q);
    if (f == "sinh") return Q::sinh(q);
    if (f == "cosh") return Q::cosh(q);
    if (f == "tanh") return Q::tanh(q);
    if (f == "asin") return Q::asin(q);
    if (f == "acos") return Q::acos(q);
    if (f == "atan") return Q::atan(q);
    if (f == "asinh") return Q::asinh(q);
    if (f == "acosh") return Q::acosh(q);
    if (f == "atanh") return Q::atanh(q);
    if (f == "tgamma") return Q::tgamma(q);
    if (f == "lgamma") return Q::lgamma(q);
    if (f == "digamma") return Q::digamma(q);
    if (f == "erf") return Q::erf(q);
    if (f == "erfc") return Q::erfc(q);
    if (f == "riemann_zeta") return Q::riemann_zeta(q);
    if (f == "hurwitz_zeta") return Q::hurwitz_zeta(q, p1);
    if (f == "polylog") return Q::polylog(p1, q);
    if (f == "lambert_w") return Q::lambert_w(n1, q);
    if (f == "expint") return Q::expint(q);
    if (f == "expint_n") return Q::expint(u1, q);
    if (f == "tgamma_upper") return Q::tgamma_upper(p1, q);
    if (f == "cyl_bessel_j") return Q::cyl_bessel_j(p1, q);
    if (f == "cyl_bessel_j_int") return Q::cyl_bessel_j(n1, q);
    if (f == "cyl_neumann") return Q::cyl_neumann(p1, q);
    if (f == "cyl_bessel_i") return Q::cyl_bessel_i(p1, q);
    if (f == "cyl_bessel_k") return Q::cyl_bessel_k(p1, q);
    if (f == "cyl_hankel_1") return Q::cyl_hankel_1(p1, q);
    if (f == "sph_bessel") return Q::sph_bessel(u1, q);
    if (f == "sph_neumann") return Q::sph_neumann(u1, q);
    if (f == "airy_ai") return Q::airy_ai(q);
    if (f == "airy_bi") return Q::airy_bi(q);
    // FMA audit (fuzz/fma_audit.py): functions without a nightly bound
    const unsigned u2 = static_cast<unsigned>(p2);
    if (f == "comp_ellint_1") return Q::comp_ellint_1(q);
    if (f == "comp_ellint_2") return Q::comp_ellint_2(q);
    if (f == "comp_ellint_3") return Q::comp_ellint_3(p1, q);
    if (f == "ellint_1") return Q::ellint_1(p1, q);
    if (f == "ellint_2") return Q::ellint_2(p1, q);
    if (f == "ellint_3") return Q::ellint_3(p1, q, p2);
    if (f == "sinint") return Q::sinint(q);
    if (f == "cosint") return Q::cosint(q);
    if (f == "fresnel_c") return Q::fresnel_c(q);
    if (f == "fresnel_s") return Q::fresnel_s(q);
    if (f == "gamma_p") return Q::gamma_p(p1, q);
    if (f == "gamma_q") return Q::gamma_q(p1, q);
    if (f == "tgamma_lower") return Q::tgamma_lower(p1, q);
    if (f == "beta") return Q::beta(p1, q);
    if (f == "polygamma") return Q::polygamma(u1, q);
    if (f == "trigamma") return Q::trigamma(q);
    if (f == "legendre") return Q::legendre(u1, q);
    if (f == "assoc_legendre") return Q::assoc_legendre(u1, u2, q);
    if (f == "laguerre") return Q::laguerre(u1, q);
    if (f == "assoc_laguerre") return Q::assoc_laguerre(u1, u2, q);
    if (f == "hermite") return Q::hermite(u1, q);
    if (f == "chebyshev_t") return Q::chebyshev_t(u1, q);
    if (f == "chebyshev_u") return Q::chebyshev_u(u1, q);
    if (f == "cyl_hankel_2") return Q::cyl_hankel_2(p1, q);
    if (f == "sph_hankel_1") return Q::sph_hankel_1(u1, q);
    const T nan = std::numeric_limits<T>::quiet_NaN();
    return {nan, nan, nan, nan};
}

template <class T>
static void run() {
    char fn[64], s1[64], s2[64], sa[64], sb[64];
    int axis;
    while (std::scanf("%63s %63s %63s %63s %63s %d", fn, s1, s2, sa, sb, &axis) == 6) {
        const T p1 = static_cast<T>(std::strtold(s1, nullptr)), p2 = static_cast<T>(std::strtold(s2, nullptr));
        const T a = static_cast<T>(std::strtold(sa, nullptr)), b = static_cast<T>(std::strtold(sb, nullptr));
        const int ax = std::abs(axis);
        const T sb_ = axis < 0 ? -b : b;  // component on the (possibly reversed) axis
        quat::quaternion<T> q(a, ax == 1 ? sb_ : T(), ax == 2 ? sb_ : T(), ax == 3 ? sb_ : T());
        if (ax == 4) {
            const T c = sb_ / T(7);
            q = quat::quaternion<T>(a, T(2) * c, T(3) * c, T(6) * c);
        }
        const auto r = eval<T>(fn, p1, p2, q);
        T comp = ax == 1 ? r.R_component_2() : ax == 2 ? r.R_component_3() : r.R_component_4();
        if (ax == 4)
            comp = (T(2) * r.R_component_2() + T(3) * r.R_component_3() + T(6) * r.R_component_4()) / T(7);
        std::printf("%La %La\n", static_cast<long double>(r.real()), static_cast<long double>(axis < 0 ? -comp : comp));
    }
}

int main(int argc, char** argv) {
    if (argc > 1 && !std::strcmp(argv[1], "long double")) run<long double>();
    else run<double>();
}
