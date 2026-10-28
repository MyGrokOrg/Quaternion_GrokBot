// Diagnostic: cost per call of std::scalbn, std::logb and the library's run-time scalbn
// (quat::detail::scalbn_run), float and double, independent calls and a dependent chain.
#include <quaternion.hpp>

#include <chrono>
#include <cstdio>
#include <vector>

namespace {
using clk = std::chrono::steady_clock;
volatile double sink;

template <class T, class F>
double ns_per_call(F f, const std::vector<T>& xs, double min_time) {
    long long n = 1 << 16;
    for (;;) {
        const auto t0 = clk::now();
        T acc = 0;
        for (long long i = 0; i < n; ++i) acc += f(xs[static_cast<std::size_t>(i) & 63], static_cast<int>(i & 7) - 3);
        const double dt = std::chrono::duration<double>(clk::now() - t0).count();
        sink = static_cast<double>(acc);
        if (dt >= min_time) return dt / static_cast<double>(n) * 1e9;
        n *= 2;
    }
}
template <class T, class F>
double ns_chain(F f, T x0, double min_time) {
    long long n = 1 << 16;
    for (;;) {
        T x = x0;
        const auto t0 = clk::now();
        for (long long i = 0; i < n; ++i) x = f(x, (i & 1) ? -5 : 5);
        const double dt = std::chrono::duration<double>(clk::now() - t0).count();
        sink = static_cast<double>(x);
        if (dt >= min_time) return dt / static_cast<double>(n) * 1e9;
        n *= 2;
    }
}
template <class T>
void run(const char* name, double mt) {
    std::vector<T> xs;
    unsigned long long s = static_cast<unsigned long long>(clk::now().time_since_epoch().count());
    for (int i = 0; i < 64; ++i) { s = s * 6364136223846793005ULL + 1442695040888963407ULL; xs.push_back(static_cast<T>(1.0 + double(s >> 40) / 16777216.0)); }
    const auto sc = [](T x, int n) { return std::scalbn(x, n); };
    const auto lb = [](T x, int) { return std::logb(x); };
    const auto fs = [](T x, int n) { return quat::detail::scalbn_run(x, n); };
    std::printf("%-7s std::scalbn  independent %7.2f ns  chain %7.2f ns\n", name, ns_per_call<T>(sc, xs, mt), ns_chain<T>(sc, T(1.25), mt));
    std::printf("%-7s std::logb    independent %7.2f ns\n", name, ns_per_call<T>(lb, xs, mt));
    std::printf("%-7s scalbn_run   independent %7.2f ns  chain %7.2f ns\n", name, ns_per_call<T>(fs, xs, mt), ns_chain<T>(fs, T(1.25), mt));
}
}  // namespace
int main() {
    run<float>("float", 0.3);
    run<double>("double", 0.3);
}
