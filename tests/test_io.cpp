#include "test_framework.hpp"
#include <format>
#include <iomanip>
#include <sstream>
using quat::quaternion;

template <class T> void stream_out() {
    std::ostringstream os;
    os << quaternion<T>(1, -2.5, 3, 0.125);
    CHECK(os.str() == "(1,-2.5,3,0.125)");
    std::ostringstream fixed;
    fixed << std::fixed << std::setprecision(2) << quaternion<T>(1, 2, 3, 4);
    CHECK(fixed.str() == "(1.00,2.00,3.00,4.00)");
    std::ostringstream w;
    w << std::setw(12) << quaternion<T>(1, 2, 3, 4) << '|';
    CHECK(w.str() == "   (1,2,3,4)|");  // width applies to the whole, as std::complex
    std::wostringstream ws;
    ws << quaternion<T>(1, 2, 3, 4);
    CHECK(ws.str() == L"(1,2,3,4)");
}
TEST(stream_output) { FOR_EACH_TYPE(stream_out); }

template <class T> void stream_in() {
    using Q = quaternion<T>;
    auto parse = [](const char* s, Q& q) { std::istringstream is(s); is >> q; return !is.fail(); };
    Q q;
    CHECK(parse("(1,2,3,4)", q) && q == Q(1, 2, 3, 4));
    CHECK(parse(" ( 1 , -2 , 3.5 , 4e1 ) ", q) && q == Q(1, -2, 3.5, 40));
    CHECK(parse("(5)", q) && q == Q(5));
    CHECK(parse("(5,6)", q) && q == Q(5, 6));
    CHECK(parse("(5,6,7)", q) && q == Q(5, 6, 7));
    CHECK(parse("7.5", q) && q == Q(7.5));
    q = Q(9, 9, 9, 9);
    CHECK(!parse("(1,2,3,4,5)", q) && q == Q(9, 9, 9, 9));
    CHECK(!parse("(1;2)", q) && q == Q(9, 9, 9, 9));
    CHECK(!parse("(1,2", q));
    CHECK(!parse("abc", q));
    // round trip
    const Q orig(T(0.1), T(-1e10), T(3.25), T(-0.0));
    std::stringstream ss;
    ss << std::setprecision(std::numeric_limits<T>::max_digits10) << orig;
    Q back; ss >> back;
    CHECK(back == orig && std::signbit(back.z()));
    std::wistringstream wis(L"(1,2,3,4)");
    Q wq; wis >> wq;
    CHECK(wq == Q(1, 2, 3, 4));
}
TEST(stream_input) { FOR_EACH_TYPE(stream_in); }

template <class T> void formatting() {
    const quaternion<T> q(1, -2.5, 3, 0.125);
    CHECK(std::format("{}", q) == "(1,-2.5,3,0.125)");
    CHECK(std::format("{:.2f}", q) == "(1.00,-2.50,3.00,0.12)");
    CHECK(std::format("{:+}", quaternion<T>(1, 2, -3, 0)) == "(+1,+2,-3,+0)");
    CHECK(std::format("{:e}", quaternion<T>(1)) == "(1.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00)");
    CHECK(std::format("{:>4}", quaternion<T>(1, 2, 3, 4)) == "(   1,   2,   3,   4)");
    CHECK(std::format(L"{}", q) == L"(1,-2.5,3,0.125)");
    CHECK(std::format("{0} {0:.1f}", quaternion<T>(0.25)) == "(0.25,0,0,0) (0.2,0.0,0.0,0.0)");
}
TEST(std_format) { FOR_EACH_TYPE(formatting); }

TEST_MAIN
