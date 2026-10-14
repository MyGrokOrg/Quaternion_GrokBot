// Minimal self-contained test harness (no external dependencies).
#pragma once
#include <cmath>
#include <cstdio>
#include <functional>
#include <limits>
#include <source_location>
#include <string>
#include <vector>
#include <quaternion/quaternion.hpp>

namespace qtest {
inline int failures = 0, checks = 0;
struct registry { std::vector<std::pair<std::string, std::function<void()>>> tests; };
inline registry& reg() { static registry r; return r; }
struct add { add(const char* n, std::function<void()> f) { reg().tests.emplace_back(n, std::move(f)); } };

inline void check(bool ok, const char* expr, std::source_location l = std::source_location::current()) {
    ++checks;
    if (!ok) { ++failures; std::fprintf(stderr, "%s:%u: CHECK failed: %s\n", l.file_name(), l.line(), expr); }
}
template <class T> T tol() {
    return std::numeric_limits<T>::epsilon() * T(256);
}
template <class T> bool near(T a, T b, T rel) {
    if (std::isnan(a) || std::isnan(b)) return std::isnan(a) && std::isnan(b);
    return std::abs(a - b) <= rel * std::max(T(1), std::max(std::abs(a), std::abs(b)));
}
template <class T> bool near(const quat::quaternion<T>& p, const quat::quaternion<T>& q, T rel) {
    return near(p.w(), q.w(), rel) && near(p.x(), q.x(), rel) && near(p.y(), q.y(), rel) && near(p.z(), q.z(), rel);
}
inline int run() {
    for (auto& [n, f] : reg().tests) {
        int before = failures;
        f();
        std::printf("[%s] %s\n", failures == before ? " OK " : "FAIL", n.c_str());
    }
    std::printf("%d checks, %d failures\n", checks, failures);
    return failures ? 1 : 0;
}
}  // namespace qtest

#define QT_CAT2(a, b) a##b
#define QT_CAT(a, b) QT_CAT2(a, b)
#define TEST(name) static void QT_CAT(test_, name)(); static qtest::add QT_CAT(reg_, name)(#name, QT_CAT(test_, name)); static void QT_CAT(test_, name)()
#define CHECK(...) qtest::check(static_cast<bool>(__VA_ARGS__), #__VA_ARGS__)
#define CHECK_NEAR(a, b, rel) qtest::check(qtest::near((a), (b), (rel)), #a " ~= " #b)
// Run a templated body for float, double and long double.
#define FOR_EACH_TYPE(fn) do { fn<float>(); fn<double>(); fn<long double>(); } while (0)
#define TEST_MAIN int main() { return qtest::run(); }
