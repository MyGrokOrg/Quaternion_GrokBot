// Evaluator for the nightly mpmath fuzz job (fuzz/fuzz.py).  Reads lines
//   fn p1 p2 a b axis
// (hex floats; axis 1..3 = i, j, k, negative = reversed) and evaluates fn on the quaternion
// a + b e_axis in the type given by argv[1] ("double" or "long double").  Prints the real part
// and the component along e_axis in hex, or "nan nan" for an unknown function.
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
    (void)p2;
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
        const auto r = eval<T>(fn, p1, p2, q);
        const T comp = ax == 1 ? r.R_component_2() : ax == 2 ? r.R_component_3() : r.R_component_4();
        std::printf("%La %La\n", static_cast<long double>(r.real()), static_cast<long double>(axis < 0 ? -comp : comp));
    }
}

int main(int argc, char** argv) {
    if (argc > 1 && !std::strcmp(argv[1], "long double")) run<long double>();
    else run<double>();
}
