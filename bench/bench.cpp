// Self-contained micro-benchmarks (no external dependencies).
//   bench [--filter SUBSTR] [--min-time SECONDS] [--json FILE]
// Each case is run on a small fixed set of inputs, repeated until --min-time has
// elapsed (default 0.2 s); the reported figure is the mean time per call.
#include <quaternion.hpp>

#include <chrono>
#include <cstdio>
#include <cstring>
#include <functional>
#include <string>
#include <vector>

namespace {
using clk = std::chrono::steady_clock;
template <class T>
using Q = quat::quaternion<T>;

volatile double sink;  // defeats dead-code elimination

struct bench_case {
    std::string group, name, type;
    std::function<double(long)> run;  // runs n calls, returns a checksum
};
std::vector<bench_case>& cases() {
    static std::vector<bench_case> c;
    return c;
}

template <class T>
const char* tname() {
    if constexpr (std::is_same_v<T, float>) return "float";
    else if constexpr (std::is_same_v<T, double>) return "double";
    else return "long double";
}

// Register f applied to each of the quaternions in `in`.
template <class T, class F>
void add(const char* group, const char* name, std::vector<Q<T>> in, F f) {
    cases().push_back({group, name, tname<T>(), [in = std::move(in), f](long n) {
                           double acc = 0;
                           const std::size_t m = in.size();
                           for (long i = 0; i < n; ++i) {
                               const auto r = f(in[static_cast<std::size_t>(i) % m]);
                               acc += static_cast<double>(quat::real(r));
                           }
                           return acc;
                       }});
}

template <class T>
std::vector<Q<T>> around(T scale) {
    std::vector<Q<T>> v;
    const T s[][4] = {{0.75, 0.5, -0.25, 0.125}, {-1.25, 0.25, 0.5, -0.75}, {0.3, -1.1, 0.4, 0.2},
                      {1.9, 0.05, -0.6, 0.35}, {-0.4, 0.8, 0.9, -0.1}, {0.15, -0.2, 0.05, 1.3}};
    for (auto& r : s) v.emplace_back(r[0] * scale, r[1] * scale, r[2] * scale, r[3] * scale);
    return v;
}
// Points a + i t (lifted to an axis) with t spread around `t0`.
template <class T>
std::vector<Q<T>> axis_pts(T a, T t0) {
    std::vector<Q<T>> v;
    for (int k = 0; k < 6; ++k) {
        const T t = t0 * (T(1) + T(k) / T(7));
        v.emplace_back(a, t * T(2) / T(3), -t / T(3), t * T(2) / T(3));
    }
    return v;
}

template <class T>
void core() {
    const auto in = around<T>(T(1));
    const Q<T> p(T(0.5), T(-0.25), T(1.5), T(0.75));
    add<T>("core", "mul", in, [p](const Q<T>& q) { return q * p; });
    add<T>("core", "div", in, [p](const Q<T>& q) { return q / p; });
    add<T>("core", "add", in, [p](const Q<T>& q) { return q + p; });
    add<T>("core", "abs", in, [](const Q<T>& q) { return Q<T>(quat::abs(q)); });
    add<T>("core", "inverse", in, [](const Q<T>& q) { return quat::inverse(q); });
    add<T>("core", "rotate", in, [p](const Q<T>& q) { return p * q * quat::conj(p); });
    add<T>("core", "scalar_div", in, [](const Q<T>& q) { return T(1.5) / q; });
    add<T>("core", "norm", in, [](const Q<T>& q) { return Q<T>(quat::norm(q)); });
    add<T>("core", "normalize", in, [](const Q<T>& q) { return quat::normalize(q); });
    add<T>("core", "proj", in, [](const Q<T>& q) { return quat::proj(q); });
    add<T>("core", "pow_int", in, [](const Q<T>& q) { return quat::pow(q, 5); });
    add<T>("core", "pow_neg_int", in, [](const Q<T>& q) { return quat::pow(q, -3); });
}
template <class T>
void elementary() {
    const auto in = around<T>(T(1));
#define E(fn) add<T>("elementary", #fn, in, [](const Q<T>& q) { return quat::fn(q); })
    E(exp); E(log); E(log10); E(sqrt); E(sin); E(cos); E(tan); E(sinh); E(cosh); E(tanh);
    E(asin); E(acos); E(atan); E(asinh); E(acosh); E(atanh);
#undef E
    add<T>("elementary", "pow", in, [](const Q<T>& q) { return quat::pow(q, Q<T>(T(0.5), T(0.25), T(0), T(0.5))); });
}
template <class T>
void special() {
    const auto in = around<T>(T(1));
    add<T>("special", "tgamma", in, [](const Q<T>& q) { return quat::tgamma(q); });
    add<T>("special", "lgamma", in, [](const Q<T>& q) { return quat::lgamma(q); });
    add<T>("special", "erf", in, [](const Q<T>& q) { return quat::erf(q); });
    add<T>("special", "erf(0.8+1.5i)", axis_pts<T>(T(0.8), T(1.5)), [](const Q<T>& q) { return quat::erf(q); });
    add<T>("special", "erfc(0.6+0.5i)", axis_pts<T>(T(0.6), T(0.5)), [](const Q<T>& q) { return quat::erfc(q); });
    add<T>("special", "erfc(2.5+i)", axis_pts<T>(T(2.5), T(1)), [](const Q<T>& q) { return quat::erfc(q); });
    add<T>("special", "erfc(0.6+5i)", axis_pts<T>(T(0.6), T(5)), [](const Q<T>& q) { return quat::erfc(q); });
    add<T>("special", "expint", in, [](const Q<T>& q) { return quat::expint(q); });
    add<T>("special", "lambert_w", in, [](const Q<T>& q) { return quat::lambert_w(q); });
    add<T>("special", "airy_ai", in, [](const Q<T>& q) { return quat::airy_ai(q); });
    add<T>("special", "gamma_upper(1.5)", in, [](const Q<T>& q) { return quat::tgamma_upper(T(1.5), q); });
    add<T>("special", "comp_ellint_1", in, [](const Q<T>& q) { return quat::comp_ellint_1(q); });
    add<T>("special", "comp_ellint_2", in, [](const Q<T>& q) { return quat::comp_ellint_2(q); });
    add<T>("special", "comp_ellint_3(0.6)", in, [](const Q<T>& q) { return quat::comp_ellint_3(T(0.6), q); });
    add<T>("special", "ellint_1(0.8)", in, [](const Q<T>& q) { return quat::ellint_1(T(0.8), q); });
    add<T>("special", "ellint_2(0.8)", in, [](const Q<T>& q) { return quat::ellint_2(T(0.8), q); });
    add<T>("special", "ellint_3(0.6,1.2)", in, [](const Q<T>& q) { return quat::ellint_3(T(0.6), q, T(1.2)); });
    // Bessel: small, large |z| and large order
    const auto b1 = around<T>(T(2)), b30 = around<T>(T(30)), b300 = around<T>(T(300)), b3k = around<T>(T(3000));
    add<T>("bessel", "J_2.5(|z|~2)", b1, [](const Q<T>& q) { return quat::cyl_bessel_j(T(2.5), q); });
    const auto b05 = around<T>(T(0.5)), b6 = around<T>(T(6)), b10 = around<T>(T(10));
    add<T>("bessel", "J_2.5(|z|~0.5)", b05, [](const Q<T>& q) { return quat::cyl_bessel_j(T(2.5), q); });
    add<T>("bessel", "J_2.5(|z|~6)", b6, [](const Q<T>& q) { return quat::cyl_bessel_j(T(2.5), q); });
    add<T>("bessel", "J_0.3(|z|~10)", b10, [](const Q<T>& q) { return quat::cyl_bessel_j(T(0.3), q); });
    add<T>("bessel", "J_3(|z|~6)", b6, [](const Q<T>& q) { return quat::cyl_bessel_j(3, q); });
    add<T>("bessel", "I_1.5(|z|~2)", b1, [](const Q<T>& q) { return quat::cyl_bessel_i(T(1.5), q); });
    add<T>("bessel", "I_0.3(|z|~6)", b6, [](const Q<T>& q) { return quat::cyl_bessel_i(T(0.3), q); });
    add<T>("bessel", "sph_j_2(|z|~2)", b1, [](const Q<T>& q) { return quat::sph_bessel(2u, q); });
    add<T>("bessel", "Y_1.5(|z|~2)", b1, [](const Q<T>& q) { return quat::cyl_neumann(T(1.5), q); });
    add<T>("bessel", "K_0.3(|z|~2)", b1, [](const Q<T>& q) { return quat::cyl_bessel_k(T(0.3), q); });
    add<T>("bessel", "Y_0.3(|z|~6)", b6, [](const Q<T>& q) { return quat::cyl_neumann(T(0.3), q); });
    add<T>("bessel", "K_0.5(|z|~2)", b1, [](const Q<T>& q) { return quat::cyl_bessel_k(T(0.5), q); });
    add<T>("bessel", "J_2.5(|z|~30)", b30, [](const Q<T>& q) { return quat::cyl_bessel_j(T(2.5), q); });
    add<T>("bessel", "J_2.5(|z|~300)", b300, [](const Q<T>& q) { return quat::cyl_bessel_j(T(2.5), q); });
    add<T>("bessel", "J_2.5(|z|~3000)", b3k, [](const Q<T>& q) { return quat::cyl_bessel_j(T(2.5), q); });
    add<T>("bessel", "J_3(|z|~300)", b300, [](const Q<T>& q) { return quat::cyl_bessel_j(3, q); });
    add<T>("bessel", "J_3(|z|~3000)", b3k, [](const Q<T>& q) { return quat::cyl_bessel_j(3, q); });
    add<T>("bessel", "Y_1.5(|z|~300)", b300, [](const Q<T>& q) { return quat::cyl_neumann(T(1.5), q); });
    add<T>("bessel", "K_0.5(|z|~300)", b300, [](const Q<T>& q) { return quat::cyl_bessel_k(T(0.5), q); });
    add<T>("bessel", "I_1.5(|z|~300)", b300, [](const Q<T>& q) { return quat::cyl_bessel_i(T(1.5), q); });
    add<T>("bessel", "J_50(|z|~30)", b30, [](const Q<T>& q) { return quat::cyl_bessel_j(T(50), q); });
    add<T>("bessel", "J_200(|z|~300)", b300, [](const Q<T>& q) { return quat::cyl_bessel_j(T(200), q); });
    add<T>("bessel", "K_100(|z|~30)", b30, [](const Q<T>& q) { return quat::cyl_bessel_k(T(100), q); });
    // zeta: small and large |Im s| on and off the critical line
    add<T>("zeta", "zeta(|s|~1)", in, [](const Q<T>& q) { return quat::riemann_zeta(q); });
    add<T>("zeta", "zeta(1/2+i 1e2)", axis_pts<T>(T(0.5), T(1e2)), [](const Q<T>& q) { return quat::riemann_zeta(q); });
    add<T>("zeta", "zeta(1/2+i 1e3)", axis_pts<T>(T(0.5), T(1e3)), [](const Q<T>& q) { return quat::riemann_zeta(q); });
    add<T>("zeta", "zeta(1/2+i 1e4)", axis_pts<T>(T(0.5), T(1e4)), [](const Q<T>& q) { return quat::riemann_zeta(q); });
    add<T>("zeta", "zeta(2+i 1e4)", axis_pts<T>(T(2), T(1e4)), [](const Q<T>& q) { return quat::riemann_zeta(q); });
    add<T>("zeta", "zeta(1/2+i 1e5)", axis_pts<T>(T(0.5), T(1e5)), [](const Q<T>& q) { return quat::riemann_zeta(q); });
    add<T>("zeta", "hurwitz(-1/2+20i,2.75)", axis_pts<T>(T(-0.5), T(20)), [](const Q<T>& q) { return quat::hurwitz_zeta(q, T(2.75)); });
    add<T>("zeta", "hurwitz(-7/4+120i,1.25)", axis_pts<T>(T(-1.75), T(120)), [](const Q<T>& q) { return quat::hurwitz_zeta(q, T(1.25)); });
    add<T>("zeta", "polylog(1/2,-3+i)", axis_pts<T>(T(-3), T(1)), [](const Q<T>& q) { return quat::polylog(T(0.5), q); });
}
}  // namespace

int main(int argc, char** argv) {
    std::string filter, json;
    double min_time = 0.2;
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--filter") && i + 1 < argc) filter = argv[++i];
        else if (!std::strcmp(argv[i], "--min-time") && i + 1 < argc) min_time = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "--json") && i + 1 < argc) json = argv[++i];
        else { std::fprintf(stderr, "usage: %s [--filter S] [--min-time SEC] [--json FILE]\n", argv[0]); return 2; }
    }
    core<float>(); core<double>(); core<long double>();
    elementary<double>(); elementary<long double>();
    special<double>(); special<long double>();
    FILE* js = json.empty() ? nullptr : std::fopen(json.c_str(), "w");
    if (js) std::fprintf(js, "{\n  \"benchmarks\": [\n");
    bool first = true;
    std::printf("%-12s %-22s %-12s %14s %12s\n", "group", "case", "type", "ns/call", "calls");
    for (auto& c : cases()) {
        const std::string full = c.group + "/" + c.name + "/" + c.type;
        if (!filter.empty() && full.find(filter) == std::string::npos) continue;
        long n = 1;
        double secs = 0;
        for (;;) {
            const auto t0 = clk::now();
            sink = c.run(n);
            secs = std::chrono::duration<double>(clk::now() - t0).count();
            if (secs >= min_time || n > (1L << 40)) break;
            n = secs > 0 ? std::max(2 * n, static_cast<long>(n * min_time / secs * 1.2)) : 2 * n;
        }
        const double ns = secs * 1e9 / static_cast<double>(n);
        std::printf("%-12s %-22s %-12s %14.1f %12ld\n", c.group.c_str(), c.name.c_str(), c.type.c_str(), ns, n);
        std::fflush(stdout);
        if (js) {
            std::fprintf(js, "%s    {\"name\": \"%s\", \"ns_per_call\": %.3f, \"calls\": %ld}", first ? "" : ",\n",
                         full.c_str(), ns, n);
            first = false;
        }
    }
    if (js) { std::fprintf(js, "\n  ]\n}\n"); std::fclose(js); }
}
