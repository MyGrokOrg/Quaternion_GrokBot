#include "test_framework.hpp"
#include <tuple>
using quat::quaternion;

// ---- compile-time checks ----
static_assert(std::is_same_v<quaternion<float>::value_type, float>);
static_assert(std::is_trivially_copyable_v<quaternion<double>>);
static_assert(sizeof(quaternion<double>) == 4 * sizeof(double));
static_assert(std::is_nothrow_copy_constructible_v<quaternion<double>>);
static_assert(std::tuple_size_v<quaternion<double>> == 4);
static_assert(std::is_same_v<std::tuple_element_t<2, quaternion<float>>, float>);
// Widening is implicit, narrowing explicit (as std::complex).
static_assert(std::is_convertible_v<quaternion<float>, quaternion<double>>);
static_assert(!std::is_convertible_v<quaternion<double>, quaternion<float>>);
static_assert(std::is_constructible_v<quaternion<float>, quaternion<double>>);
static_assert(std::is_convertible_v<std::complex<double>, quaternion<double>>);
static_assert(std::is_same_v<decltype(quaternion(1.0f, 2.0f, 3.0f, 4.0f)), quaternion<float>>);
static_assert(std::is_same_v<decltype(quaternion(std::complex<double>{})), quaternion<double>>);

constexpr quaternion<double> I(0, 1, 0, 0), J(0, 0, 1, 0), K(0, 0, 0, 1);
static_assert(I * I == -1.0 && J * J == -1.0 && K * K == -1.0);
static_assert(I * J * K == -1.0);
static_assert(I * J == K && J * K == I && K * I == J);
static_assert(J * I == -K && K * J == -I && I * K == -J);
static_assert(quaternion<double>(1, 1, 1, 1) / quaternion<double>(1, -1, 1, -1) == quaternion<double>(0, 1, 0, 0));
static_assert(conj(quaternion<double>(1, 2, 3, 4)) == quaternion<double>(1, -2, -3, -4));
static_assert(norm(quaternion<double>(1, 2, 3, 4)) == 30.0);
static_assert(pow(I, 2) == -1.0 && pow(I, 4) == 1.0 && pow(I, -1) == -I && pow(I, 0) == 1.0);
static_assert(dot(quaternion<double>(1, 2, 3, 4), quaternion<double>(5, 6, 7, 8)) == 70.0);
static_assert(cross(I, J) == K);
static_assert(inverse(quaternion<double>(2)) == 0.5);
static_assert([] { quaternion<double> q(1, 2, 3, 4); q.get<3>() = 9; q.x(7); return q; }() ==
              quaternion<double>(1, 7, 3, 9));
static_assert(quat::get<1>(quaternion<double>(1, 2, 3, 4)) == 2.0);

template <class T> void ctor_access() {
    using Q = quaternion<T>;
    Q d;
    CHECK(d.real() == T(0) && d.R_component_2() == 0 && d.R_component_3() == 0 && d.R_component_4() == 0);
    Q q(1, 2, 3, 4);
    CHECK(q.real() == 1 && q.R_component_1() == 1 && q.R_component_2() == 2 && q.R_component_3() == 3 && q.R_component_4() == 4);
    CHECK(q.w() == 1 && q.x() == 2 && q.y() == 3 && q.z() == 4);
    CHECK(real(q) == 1);
    CHECK(unreal(q) == Q(0, 2, 3, 4) && q.unreal() == Q(0, 2, 3, 4));
    CHECK(q.C_component_1() == std::complex<T>(1, 2) && q.C_component_2() == std::complex<T>(3, 4));
    Q fromc(std::complex<T>(1, 2), std::complex<T>(3, 4));
    CHECK(fromc == q);
    Q one(std::complex<T>(5, 6));
    CHECK(one == Q(5, 6, 0, 0));
    q.real(9); CHECK(q.w() == 9);
    q.w(1); q.x(-1); q.y(-2); q.z(-3);
    CHECK(q == Q(1, -1, -2, -3));
    Q a; a = T(3); CHECK(a == Q(3));
    a = std::complex<T>(1, 2); CHECK(a == Q(1, 2));
    a = quaternion<float>(1, 2, 3, 4); CHECK(a == Q(1, 2, 3, 4));
    quaternion<long double> wide = quaternion<float>(1.5f, 2, 3, 4);
    CHECK(wide == quaternion<long double>(1.5L, 2, 3, 4));
    quaternion<float> narrow(quaternion<double>(0.25, 1, 2, 3));
    CHECK(narrow == quaternion<float>(0.25f, 1, 2, 3));
}
TEST(constructors_and_accessors) { FOR_EACH_TYPE(ctor_access); }

template <class T> void tuple_like() {
    quaternion<T> q(1, 2, 3, 4);
    auto [a, b, c, d] = q;
    CHECK(a == 1 && b == 2 && c == 3 && d == 4);
    auto& [ra, rb, rc, rd] = q;
    rd = 8; (void)ra; (void)rb; (void)rc;
    CHECK(q.z() == 8);
    using std::get;  // ADL finds quat::get
    CHECK(get<0>(q) == 1 && get<3>(q) == 8);
    static_assert(std::is_same_v<decltype(get<0>(std::move(q))), T&&>);
}
TEST(tuple_protocol) { FOR_EACH_TYPE(tuple_like); }

template <class T> void arithmetic() {
    using Q = quaternion<T>;
    using C = std::complex<T>;
    const Q p(1, 2, 3, 4), q(5, -6, 7, -8);
    CHECK(+p == p);
    CHECK(-p == Q(-1, -2, -3, -4));
    CHECK(p + q == Q(6, -4, 10, -4));
    CHECK(p - q == Q(-4, 8, -4, 12));
    // Hamilton product computed by hand.
    CHECK(p * q == Q(1*5 - 2*-6 - 3*7 - 4*-8, 1*-6 + 2*5 + 3*-8 - 4*7, 1*7 - 2*-8 + 3*5 + 4*-6, 1*-8 + 2*7 - 3*-6 + 4*5));
    CHECK(p * q != q * p);  // non-commutative
    CHECK(norm(p * q) == norm(p) * norm(q));  // multiplicative norm
    CHECK(conj(p * q) == conj(q) * conj(p));
    CHECK_NEAR((p * q) / q, p, qtest::tol<T>());
    CHECK_NEAR(p / q * q, p, qtest::tol<T>());
    CHECK_NEAR(p * inverse(p), Q(1), qtest::tol<T>());
    CHECK_NEAR(inverse(p) * p, Q(1), qtest::tol<T>());
    CHECK((p * q) * Q(0, 1, 1, 0) == p * (q * Q(0, 1, 1, 0)));  // associativity (exact for small ints)
    // scalar
    CHECK(p + T(1) == Q(2, 2, 3, 4) && T(1) + p == Q(2, 2, 3, 4));
    CHECK(p - T(1) == Q(0, 2, 3, 4) && T(1) - p == Q(0, -2, -3, -4));
    CHECK(p * T(2) == Q(2, 4, 6, 8) && T(2) * p == Q(2, 4, 6, 8));
    CHECK(p / T(2) == Q(T(0.5), 1, T(1.5), 2));
    CHECK_NEAR(T(2) / p, T(2) * inverse(p), qtest::tol<T>());
    CHECK(p * 2 == Q(2, 4, 6, 8));  // int converts through hidden friends
    // complex
    const C z(1, 1);
    CHECK(p + z == Q(2, 3, 3, 4) && z + p == Q(2, 3, 3, 4));
    CHECK(p - z == Q(0, 1, 3, 4) && z - p == Q(0, -1, -3, -4));
    CHECK(p * z == p * Q(z) && z * p == Q(z) * p && p * z != z * p);
    CHECK_NEAR(p / z, p * inverse(Q(z)), qtest::tol<T>());
    CHECK_NEAR(z / p, Q(z) * inverse(p), qtest::tol<T>());
    // complex numbers embed faithfully
    CHECK(Q(C(1, 2)) * Q(C(3, 4)) == Q(C(1, 2) * C(3, 4)));
    // compound
    Q r = p; r += q; CHECK(r == p + q);
    r = p; r -= q; CHECK(r == p - q);
    r = p; r *= q; CHECK(r == p * q);
    r = p; r /= q; CHECK(r == p / q);
    r = p; r += T(1); r -= T(2); r *= T(4); r /= T(2); CHECK(r == Q(0, 4, 6, 8));
    r = p; r += z; r -= z; r *= z; r /= z; CHECK_NEAR(r, p, qtest::tol<T>());
    r = p; r *= r; CHECK(r == p * p);  // self-aliasing
    r = p; r += quaternion<float>(1, 1, 1, 1); CHECK(r == Q(2, 3, 4, 5));
}
TEST(arithmetic) { FOR_EACH_TYPE(arithmetic); }

template <class T> void equality() {
    using Q = quaternion<T>;
    CHECK(Q(1) == T(1) && T(1) == Q(1) && Q(1, 1) != T(1));
    CHECK(Q(1, 2) == std::complex<T>(1, 2) && std::complex<T>(1, 2) == Q(1, 2));
    CHECK(Q(1, 2, 3) != std::complex<T>(1, 2));
    CHECK(Q(0.0) == Q(-0.0));
    const T nan = std::numeric_limits<T>::quiet_NaN();
    CHECK(Q(nan) != Q(nan));
}
TEST(equality) { FOR_EACH_TYPE(equality); }

template <class T> void helpers() {
    using Q = quaternion<T>;
    const Q p(1, 2, 3, 4), q(-2, 1, 0, 5);
    CHECK(dot(p, q) == T(-2 + 2 + 0 + 20));
    CHECK(cross(p, q) == Q(0, 3*5 - 4*0, 4*1 - 2*5, 2*0 - 3*1));
    // For pure quaternions u*v = -dot + cross.
    const Q u = p.unreal(), v = q.unreal();
    CHECK(u * v == -dot(u, v) + cross(u, v));
    CHECK_NEAR(abs(normalize(p)), T(1), qtest::tol<T>());
    CHECK(normalize(Q()) == Q());
}
TEST(helpers) { FOR_EACH_TYPE(helpers); }

template <class T> void literals() {
    using namespace quat::literals;
    CHECK(2.0_i == quaternion<double>(0, 2, 0, 0) && 3_j == quaternion<double>(0, 0, 3, 0) && 4.5_k == quaternion<double>(0, 0, 0, 4.5));
    CHECK(1.0_if == quaternion<float>(0, 1) && 1_jf == quaternion<float>(0, 0, 1) && 1.0_kf == quaternion<float>(0, 0, 0, 1));
    CHECK(1.0_il == quaternion<long double>(0, 1) && 1_jl == quaternion<long double>(0, 0, 1) && 2_kl == quaternion<long double>(0, 0, 0, 2));
    CHECK(1.0 + 2.0_i + 3.0_j + 4.0_k == quaternion<double>(1, 2, 3, 4));
    CHECK(1_i * 1_j == 1_k);
    static_assert(1_i * 1_i == -1.0);
}
TEST(literals) { literals<double>(); }

TEST_MAIN
