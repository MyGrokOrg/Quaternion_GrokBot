#include <quaternion.hpp>
#include <iostream>
#include <numbers>
#include <print>

int main() {
    using namespace quat::literals;
    using Q = quat::quaternion<double>;

    Q p = 1.0 + 2.0_i + 3.0_j + 4.0_k;
    Q q(0.5, -1.0, 0.25, 2.0);

    std::cout << "p       = " << p << '\n';
    std::println("p * q   = {:.3f}", p * q);
    std::println("q * p   = {:.3f}  (non-commutative)", q * p);
    std::println("p / q   = {:.3f}", p / q);
    std::println("|p|     = {:.6f}, arg(p) = {:.6f}", abs(p), arg(p));
    std::println("exp(p)  = {:.4f}", exp(p));
    std::println("sqrt(p) = {:.4f}", sqrt(p));
    std::println("i*j*k   = {}", 1_i * 1_j * 1_k);

    // Rotate the vector (1, 0, 0) by 90 degrees about the z axis.
    Q r = quat::polar(1.0, std::numbers::pi / 4, 1.0_k);
    Q v = r * 1.0_i * conj(r);
    std::println("rotated = {:.3f}", v);

    // Special functions (axis reduction of the complex functions).
    std::println("tgamma(p) = {:.6g}", quat::tgamma(p));
    std::println("erf(p)    = {:.6g}", quat::erf(p));
    std::println("zeta(p)   = {:.6g}", quat::riemann_zeta(p));
    std::println("J_2(p)    = {:.6g}", quat::cyl_bessel_j(2, p));
    std::println("sin(asin(p)) = {:.6f}", sin(quat::asin(p)));

    auto [w, x, y, z] = p;  // structured bindings
    std::println("components sum = {}", w + x + y + z);
}
