// Evaluator for fuzz/log_unit_check.py (log, log10 and pow next to |q| = 1).
// Build: g++ -std=c++23 -O2 -Iinclude fuzz/log_unit_eval.cpp -o log_unit_eval
//   "log":  stdin a b c d (hex doubles); stdout Re log, Re log10 in double, then in long double,
//           then in float (inputs rounded to float)
//   "pow":  stdin a b c d p;             stdout pow(q, p) in double (4), then long double (4)
#include <quaternion/quaternion.hpp>
#include <cstdio>
#include <cstring>
int main(int argc, char** argv) {
    const bool pw = argc > 1 && std::strcmp(argv[1], "pow") == 0;
    double a, b, c, d, p = 0;
    while (std::scanf("%la %la %la %la", &a, &b, &c, &d) == 4 && (!pw || std::scanf("%la", &p) == 1)) {
        const quat::quaternion<double> q(a, b, c, d);
        const quat::quaternion<long double> Q(a, b, c, d);
        if (pw) {
            const auto r = quat::pow(q, p);
            const auto R = quat::pow(Q, static_cast<long double>(p));
            std::printf("%a %a %a %a %La %La %La %La\n", r.w(), r.x(), r.y(), r.z(), R.w(), R.x(), R.y(), R.z());
        } else {
            const quat::quaternion<float> f(static_cast<float>(a), static_cast<float>(b), static_cast<float>(c), static_cast<float>(d));
            std::printf("%a %a %La %La %a %a\n", quat::log(q).w(), quat::log10(q).w(), quat::log(Q).w(), quat::log10(Q).w(),
                        static_cast<double>(quat::log(f).w()), static_cast<double>(quat::log10(f).w()));
        }
    }
}
