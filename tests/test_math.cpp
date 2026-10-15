#include "test_framework.hpp"
#include <numbers>
using quat::quaternion;

// Any quaternion with a fixed unit axis behaves like a complex number; check
// functions against std::complex by rotating the i axis onto u.
template <class T, class QF, class CF>
void against_complex(QF qf, CF cf, T a, T b, T rel) {
    const T ux = T(2) / T(3), uy = T(-1) / T(3), uz = T(2) / T(3);  // unit vector
    const quaternion<T> q(a, b * ux, b * uy, b * uz);
    const std::complex<T> z = cf(std::complex<T>(a, b));
    const quaternion<T> expected(z.real(), z.imag() * ux, z.imag() * uy, z.imag() * uz);
    CHECK_NEAR(qf(q), expected, rel);
    // On the i axis results coincide with std::complex.
    CHECK_NEAR(qf(quaternion<T>(std::complex<T>(a, b))), quaternion<T>(z), rel);
}

template <class T> void transcendental() {
    using Q = quaternion<T>;
    const T rel = std::numeric_limits<T>::epsilon() * T(1024);
    const long double samples[][2] = {{0.5, 0.25}, {-1.25, 2}, {2, -0.75}, {0.1, 3}, {-0.3, -0.4}, {1.5, 1e-3}};
    for (auto& s : samples) {
        const T a = static_cast<T>(s[0]), b = static_cast<T>(s[1]);
#define CMP(fn) against_complex<T>([](const Q& q) { return quat::fn(q); }, [](const std::complex<T>& z) { return std::fn(z); }, a, b, rel)
        CMP(exp); CMP(log); CMP(log10); CMP(sqrt); CMP(sin); CMP(cos); CMP(tan); CMP(sinh); CMP(cosh); CMP(tanh);
#undef CMP
        against_complex<T>([](const Q& q) { return quat::pow(q, T(1.5)); }, [](const std::complex<T>& z) { return std::pow(z, T(1.5)); }, a, b, rel);
        against_complex<T>([](const Q& q) { return quat::pow(T(2.5), q); }, [](const std::complex<T>& z) { return std::pow(T(2.5), z); }, a, b, rel);
        against_complex<T>([](const Q& q) { return quat::pow(q, 3); }, [](const std::complex<T>& z) { return z * z * z; }, a, b, rel);
        CHECK_NEAR(abs(Q(std::complex<T>(a, b))), std::abs(std::complex<T>(a, b)), rel);
        CHECK_NEAR(arg(Q(std::complex<T>(a, std::abs(b)))), std::arg(std::complex<T>(a, std::abs(b))), rel);
    }
    const Q q(T(0.3), T(-0.2), T(0.5), T(0.1)), p(T(-0.4), T(0.7), T(0.2), T(-0.6));
    CHECK_NEAR(log(exp(q)), q, rel);
    CHECK_NEAR(exp(log(p)), p, rel);
    CHECK_NEAR(sqrt(p) * sqrt(p), p, rel);
    CHECK(sqrt(p).real() >= 0);
    CHECK_NEAR(sin(q) * sin(q) + cos(q) * cos(q), Q(1), rel);
    CHECK_NEAR(cosh(q) * cosh(q) - sinh(q) * sinh(q), Q(1), rel);
    CHECK_NEAR(tan(q), sin(q) / cos(q), rel);
    CHECK_NEAR(tanh(q), sinh(q) / cosh(q), rel);
    CHECK_NEAR(pow(p, 5), p * p * p * p * p, rel);
    CHECK_NEAR(pow(p, -2), inverse(p * p), rel);
    CHECK_NEAR(pow(p, T(2)), p * p, rel);
    CHECK_NEAR(pow(p, T(0.5)), sqrt(p), rel);
    CHECK_NEAR(pow(p, Q(2)), p * p, rel);
    CHECK_NEAR(pow(p, q), exp(q * log(p)), rel);
    CHECK_NEAR(pow(p, std::complex<T>(2)), p * p, rel);
    CHECK_NEAR(pow(std::complex<T>(1, 1), Q(2)), Q(std::complex<T>(1, 1) * std::complex<T>(1, 1)), rel);
    CHECK_NEAR(pow(T(2), Q(3)), Q(8), rel);
    CHECK(pow(Q(), T(2)) == Q() && pow(Q(), T(0)) == Q(1) && pow(Q(), Q()) == Q(1));
    // exp of a pure quaternion is a unit quaternion (rotation).
    CHECK_NEAR(abs(exp(q.unreal())), T(1), rel);
    // abs / arg / norm
    CHECK_NEAR(abs(Q(1, 2, 3, 4)), std::sqrt(T(30)), rel);
    CHECK(abs(Q(-3)) == 3);
    const T big = std::numeric_limits<T>::max() / 4;
    CHECK(std::isfinite(abs(Q(big, big, big, big))));  // no overflow
    CHECK_NEAR(arg(Q(1, 1, 1, 1)), std::atan2(std::sqrt(T(3)), T(1)), rel);
    CHECK_NEAR(arg(Q(-1)), quat::detail::pi_v<T>, rel);
    CHECK(arg(Q(1)) == 0);
    // negative reals: log and sqrt pick the i axis like std::complex
    CHECK_NEAR(log(Q(-1)), Q(0, quat::detail::pi_v<T>), rel);
    CHECK_NEAR(sqrt(Q(-4)), Q(0, 2), rel);
    CHECK(sqrt(Q(9)) == Q(3));
    CHECK_NEAR(sqrt(Q(-4, 0, 3, 0)), Q(std::sqrt(T(0.5)), 0, T(3) / (2 * std::sqrt(T(0.5))), 0), rel);
}
TEST(transcendental) { FOR_EACH_TYPE(transcendental); }

// log|q| next to |q| = 1 is relatively accurate (computed from |q|^2 - 1 in double length);
// before, log(hypot) gave only absolute accuracy there (e.g. 0 for the first point).
// References: mpmath at 400 bits for the exact (binary64) inputs.
template <class T> void log_near_unit() {
    using Q = quaternion<T>;
    if constexpr (std::numeric_limits<T>::digits >= 53) {
        struct pt { double a, b, c, d; long double ln, lg; };
        const pt P[] = {
            {0x1.fffffffffffc0p-1, 0x1p-23, 0, 0, 2.524354896707237777317531e-29L, 1.0963134020054066135346e-29L},
            {0x1.ff5fdbdeb6fd2p-1, 0x1.94d38707c8210p-5, -0x1.52120d939d01ep-16, -0x1.c2f75f94187f7p-21,
             -1.245280173448013895465821e-17L, -5.408183077519967678167502e-18L},
            {0x1.3333333333333p-1, 0x1.999999999999ap-1, 0, 0, 2.220446049250313154802973e-17L, 9.643274665532871428387039e-18L},
            {0x1.ffffffc476116p-1, 0x1.3f37d02c75800p-29, 0, 0, -6.93122827933508981306832e-9L, -3.010194194527000466631502e-9L},
            {-0x1.50082a256e6a2p-4, 0x1.04f4c6405b7bfp-3, 0x1.7e875f79dce02p-1, -0x1.b04b410deb98ep-1,
             0.1288874979564298891418119L, 0.05597512914879414740825717L},
        };
        const T tol = std::numeric_limits<T>::epsilon() * T(4);
        for (const pt& p : P) {
            const Q q(T(p.a), T(p.b), T(p.c), T(p.d));
            const T ln = static_cast<T>(p.ln), lg = static_cast<T>(p.lg);
            CHECK(std::abs(quat::log(q).w() - ln) <= tol * std::abs(ln));
            CHECK(std::abs(quat::log10(q).w() - lg) <= tol * std::abs(lg));
            // the complex wrapper path (Annex G kernel) agrees on the i axis
            if (p.c == 0 && p.d == 0) CHECK(std::abs(quat::log(Q(T(p.a), T(p.b))).w() - ln) <= tol * std::abs(ln));
        }
    }
    // |q| = 1 exactly: log|q| = 0 exactly
    CHECK(quat::log(Q(T(0.5), T(0.5), T(0.5), T(-0.5))).w() == T(0));
}
TEST(log_near_unit_circle) { FOR_EACH_TYPE(log_near_unit); }

template <class T> void polar_proj() {
    using Q = quaternion<T>;
    const T rel = std::numeric_limits<T>::epsilon() * T(64);
    CHECK_NEAR(quat::polar(T(2), T(0.5)), Q(std::polar(T(2), T(0.5))), rel);
    CHECK(quat::polar(T(3)) == Q(3));
    const Q axis(0, 0, 3, 4);
    const Q p = quat::polar(T(2), T(0.7), axis);
    CHECK_NEAR(abs(p), T(2), rel);
    CHECK_NEAR(arg(p), T(0.7), rel);
    CHECK_NEAR(normalize(p.unreal()), Q(0, 0, T(3) / 5, T(4) / 5), rel);
    CHECK_NEAR(quat::polar(abs(Q(1, 2, 3, 4)), arg(Q(1, 2, 3, 4)), Q(1, 2, 3, 4)), Q(1, 2, 3, 4), rel);
    const T inf = std::numeric_limits<T>::infinity();
    CHECK(proj(Q(1, 2, 3, 4)) == Q(1, 2, 3, 4));
    const Q pr = proj(Q(1, -inf, 3, -4));
    CHECK(pr.w() == inf && pr.x() == 0 && std::signbit(pr.x()) && !std::signbit(pr.y()) && std::signbit(pr.z()));
}
TEST(polar_and_proj) { FOR_EACH_TYPE(polar_proj); }

TEST(rotation_example) {
    // Rotating (1,0,0) by 90 degrees about z gives (0,1,0).
    using Q = quaternion<double>;
    const Q r = quat::polar(1.0, quat::detail::pi_v<double> / 4, Q(0, 0, 0, 1));
    const Q v = r * Q(0, 1, 0, 0) * conj(r);
    CHECK_NEAR(v, Q(0, 0, 1, 0), 1e-15);
}

TEST_MAIN
