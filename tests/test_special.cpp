#include "test_framework.hpp"
#include <quaternion.hpp>
#include <cmath>
#include <numbers>
#include <string_view>
using quat::quaternion;

#if defined(__cpp_lib_math_special_functions) || (defined(__STDCPP_MATH_SPEC_FUNCS__) && __STDCPP_MATH_SPEC_FUNCS__ >= 201003L)
#define HAVE_STD_SPECIAL_MATH 1
#else
#define HAVE_STD_SPECIAL_MATH 0
#endif

namespace {

// Relative tolerance: kernels run in at least double precision.
template <class T> T rtol(T k = T(2000)) {
    using W = quat::detail::sf::work_t<T>;
    return std::max(T(256) * std::numeric_limits<T>::epsilon(),
                    static_cast<T>(k * std::numeric_limits<W>::epsilon()));
}
template <class T> bool close(const quaternion<T>& got, const quaternion<T>& want, T tol, T scale = T()) {
    scale = std::max({scale, abs(want), T(std::numeric_limits<T>::min() * 1e10)});
    return abs(got - want) <= tol * scale;
}
#define CHECK_CLOSE_S(got, want, tol, scale) qtest::check(close((got), (want), (tol), (scale)), #got " ~= " #want)
#define CHECK_CLOSE(got, want, tol) qtest::check(close((got), (want), (tol)), #got " ~= " #want)

// A unit axis not aligned with i, j or k.
template <class T> quaternion<T> axis() { return {T(0), T(2) / 3, T(-1) / 3, T(2) / 3}; }
template <class T> quaternion<T> on_axis(T a, T b) { return T(b) * axis<T>() + a; }

struct ref_value { std::string_view fn; int n; long double re, im, wre, wim; };
// Reference values computed with mpmath at 30 digits (inputs are exact binary fractions).
constexpr ref_value refs[] = {
    {"tgamma", 0, 0.75L, 0.5L, 0.8349299659737468481748878L, -0.4063818800581324332638504L},
    {"tgamma", 0, 2.5L, 1.25L, 0.5387173208874147138683515L, 0.7516029084919344805159757L},
    {"tgamma", 0, -1.75L, 0.5L, 0.8106922465006138165545441L, -0.2766688174174664664363188L},
    {"tgamma", 0, 0.25L, 3.0L, 0.01705032393424411927268868L, -0.001596877420381335891039394L},
    {"tgamma", 0, -3.5L, 2.0L, -0.001561837432876754544743225L, 0.0004611942720843740308950255L},
    {"tgamma", 0, 5.0L, 0.125L, 23.53513590594841585977701L, 4.484323991768101778126279L},
    {"tgamma", 0, 1.5L, 0.0L, 0.8862269254527580136490837L, 0.0L},
    {"tgamma", 0, -0.5L, 6.5L, -0.000008970556558623174242112469L, -0.00001093956424338666283329023L},
    {"tgamma", 0, 10.0L, 4.0L, -152693.7870394386393140697L, 48182.09543641663191095898L},
    {"tgamma", 0, 0.0L, 1.0L, -0.1549498283018106851249551L, -0.4980156681183560427136911L},
    {"lgamma", 0, 0.75L, 0.5L, -0.07410253140811996089610494L, -0.4529718950149241177499544L},
    {"lgamma", 0, 2.5L, 1.25L, -0.07825481438511577458210629L, 0.9489117675513034962886646L},
    {"lgamma", 0, -1.75L, 0.5L, -0.1547814928586282600266802L, -6.612066054514497875581948L},
    {"lgamma", 0, 0.25L, 3.0L, -4.067219409137411985568708L, -0.09338431339316938304969317L},
    {"lgamma", 0, -3.5L, 2.0L, -6.420091394575657853388377L, -9.711907658196487230468898L},
    {"lgamma", 0, 5.0L, 0.125L, 3.176324962666451453680834L, 0.1882805870530784209776885L},
    {"lgamma", 0, 1.5L, 0.0L, -0.1207822376352452223455184L, 0.0L},
    {"lgamma", 0, -0.5L, 6.5L, -11.16598962892734085685618L, 4.025565106169746721162491L},
    {"lgamma", 0, 10.0L, 4.0L, 11.98364941340949788951092L, 9.119119418031396853718077L},
    {"lgamma", 0, 0.0L, 1.0L, -0.6509231993018563388852168L, -1.872436647262429817118853L},
    {"digamma", 0, 0.75L, 0.5L, -0.6094013558546486782665547L, 0.9818526342169299431900655L},
    {"digamma", 0, 2.5L, 1.25L, 0.8614311824530846459049656L, 0.5520531988482499536811979L},
    {"digamma", 0, -1.75L, 0.5L, 0.5709625026153088809586406L, 2.914356413025859101041191L},
    {"digamma", 0, 0.25L, 3.0L, 1.097449149522477930457005L, 1.654730547313617386820997L},
    {"digamma", 0, -3.5L, 2.0L, 1.49912081258195930842447L, 2.679572480614529243610329L},
    {"digamma", 0, 5.0L, 0.125L, 1.506498695227298951819709L, 0.02765839737927733765803538L},
    {"digamma", 0, 1.5L, 0.0L, 0.03648997397857652055902367L, 0.0L},
    {"digamma", 0, -0.5L, 6.5L, 1.882576552336166258033434L, 1.72373750326548484876763L},
    {"digamma", 0, 10.0L, 4.0L, 2.333171465077168286831019L, 0.3982425817126880381473303L},
    {"digamma", 0, 0.0L, 1.0L, 0.09465032062247697727187848L, 2.076674047468581174134051L},
    {"erf", 0, 0.75L, 0.5L, 0.84155181003633210089521L, 0.316301010060480927952821L},
    {"erf", 0, 2.5L, 1.25L, 0.9983452221740688511249073L, 0.0006757278098157495351560376L},
    {"erf", 0, -1.75L, 0.5L, -1.006617597425512493740281L, 0.01537392950966603028869575L},
    {"erf", 0, 0.25L, 3.0L, 1500.805941142063371657009L, 254.3455295566197758040481L},
    {"erf", 0, -3.5L, 2.0L, -1.000012585981940830651822L, 0.00003369461013273328798794915L},
    {"erf", 0, 5.0L, 0.125L, 0.9999999999995435288730532L, 0.000000000001493012132820743744128311L},
    {"erf", 0, 1.5L, 0.0L, 0.9661051464753107270669763L, 0.0L},
    {"erf", 0, -0.5L, 6.5L, -20975014571272012.88278074L, 150884100455778401.6660491L},
    {"erf", 0, 10.0L, 4.0L, 1.0L, -0.00000000000000000000000000000000000001665104622779648303732417L},
    {"erf", 0, 0.0L, 1.0L, 0.0L, 1.650425758797542876025338L},
    {"erfc", 0, 0.75L, 0.5L, 0.15844818996366789910479L, -0.316301010060480927952821L},
    {"erfc", 0, 2.5L, 1.25L, 0.001654777825931148875092714L, -0.0006757278098157495351560376L},
    {"erfc", 0, -1.75L, 0.5L, 2.006617597425512493740281L, -0.01537392950966603028869575L},
    {"erfc", 0, 0.25L, 3.0L, -1499.805941142063371657009L, -254.3455295566197758040481L},
    {"erfc", 0, -3.5L, 2.0L, 2.000012585981940830651822L, -0.00003369461013273328798794915L},
    {"erfc", 0, 5.0L, 0.125L, 0.0000000000004564711269468142355342469L, -0.000000000001493012132820743744128311L},
    {"erfc", 0, 1.5L, 0.0L, 0.03389485352468927293302374L, 0.0L},
    {"erfc", 0, -0.5L, 6.5L, 20975014571272013.88278074L, -150884100455778401.6660491L},
    {"erfc", 0, 10.0L, 4.0L, 0.000000000000000000000000000000000000004553981839356313957621737L, 0.00000000000000000000000000000000000001665104622779648303732417L},
    {"erfc", 0, 0.0L, 1.0L, 1.0L, -1.650425758797542876025338L},
    {"zeta", 0, 0.75L, 0.5L, -0.2401410529536000482524055L, -1.562360750614463283396273L},
    {"zeta", 0, 2.5L, 1.25L, 1.076673100473679534374074L, -0.256671578461021853201539L},
    {"zeta", 0, -1.75L, 0.5L, 0.0009822619081030396406314443L, -0.02258464972051546422193441L},
    {"zeta", 0, 0.25L, 3.0L, 0.4852981185578533691235869L, -0.0589857558159271582735397L},
    {"zeta", 0, -3.5L, 2.0L, -0.003560979964919072343271255L, 0.04262253731477640726730476L},
    {"zeta", 0, 5.0L, 0.125L, 1.036742591669069900836825L, -0.003564593008265411617677386L},
    {"zeta", 0, 1.5L, 0.0L, 2.612375348685488343348568L, 0.0L},
    {"zeta", 0, -0.5L, 6.5L, 0.8333867886259631261550898L, 0.5977286648389660727867855L},
    {"zeta", 0, 10.0L, 4.0L, 0.9990846953448826632346191L, -0.0003355381072179333459854128L},
    {"zeta", 0, 0.0L, 1.0L, 0.003300223685324102874217114L, -0.4181554491413216766892742L},
    {"zeta", 0, 0.5L, 14.25L, -0.01021320888723080752252596L, 0.09169808962101764008146596L},
    {"zeta", 0, -2.5L, 10.0L, 4.263590288889194478036571L, 1.459816619917562812899585L},
    {"zeta", 0, 3.0L, 30.0L, 0.9365853681541057681990584L, -0.1359171988090308639946303L},
    {"besselj", 0, 0.75L, 0.5L, 0.9148274060428001696858415L, -0.1800909890946812142803793L},
    {"besselj", 0, 2.5L, 1.25L, -0.266264986386348745687975L, -0.7315078461281092140728955L},
    {"besselj", 0, -7.5L, 3.0L, 2.334417964838417708551494L, 1.613924213217510851458155L},
    {"besselj", 0, 15.0L, 0.0L, -0.01422447282678077323386427L, 0.0L},
    {"besselj", 0, 0.0L, 8.0L, 427.5641157218047851773968L, 0.0L},
    {"besselj", 0, 30.5L, 2.0L, -0.05586224571887695858455363L, 0.5209075034655020880583911L},
    {"besselj", 1, 0.75L, 0.5L, 0.3822856894133987563527996L, 0.2054098665433370046773922L},
    {"besselj", 1, 2.5L, 1.25L, 0.7719997131780264252019543L, -0.3900708878111239448611259L},
    {"besselj", 1, -7.5L, 3.0L, -1.719413633448867968239744L, 2.176253984798589612961258L},
    {"besselj", 1, 15.0L, 0.0L, 0.2051040386135227611471374L, 0.0L},
    {"besselj", 1, 0.0L, 8.0L, 0.0L, 399.8731367825600982190831L},
    {"besselj", 1, 30.5L, 2.0L, -0.5405866907494848080633512L, -0.04404733505684711895234636L},
    {"besselj", 3, 0.75L, 0.5L, -0.002179498124517558086889013L, 0.01480595450779401032791866L},
    {"besselj", 3, 2.5L, 1.25L, 0.1917550434001965662991532L, 0.2839096210612467423860114L},
    {"besselj", 3, -7.5L, 3.0L, 2.1592264886331963991523L, -0.9570725841605407381277808L},
    {"besselj", 3, 15.0L, 0.0L, -0.1940182578201226345550976L, 0.0L},
    {"besselj", 3, 0.0L, 8.0L, 0.0L, -236.0752210194777179077701L},
    {"besselj", 3, 30.5L, 2.0L, 0.5387823098246839004361358L, -0.02422371667289032266227796L},
    {"besselj", -2, 0.75L, 0.5L, 0.04374293323527384779919438L, 0.08880374035819721509673462L},
    {"besselj", -2, 2.5L, 1.25L, 0.6355221187207259954616655L, 0.2348225697120214332971495L},
    {"besselj", -2, -7.5L, 3.0L, -1.73903515394915246408615L, -1.956105484808095317128353L},
    {"besselj", -2, 15.0L, 0.0L, 0.04157167797525047472014926L, 0.0L},
    {"besselj", -2, 0.0L, 8.0L, -327.595831526164760622626L, 0.0L},
    {"besselj", -2, 30.5L, 2.0L, 0.0203771213132618001497944L, -0.5214689599016484461251182L},
    {"besselj", 12, 0.75L, 0.5L, 0.0000000000001058103912632777393016346L, 0.0000000000001002492453438937244667145L},
    {"besselj", 12, 2.5L, 1.25L, 0.00000007067041063162697609769175L, -0.00000007890311883640074311599296L},
    {"besselj", 12, -7.5L, 3.0L, -0.01396108868625836673112269L, 0.007594863981089882021868499L},
    {"besselj", 12, 15.0L, 0.0L, 0.236665844054768055905717L, 0.0L},
    {"besselj", 12, 0.0L, 8.0L, 0.1141962709681452775014538L, 0.0L},
    {"besselj", 12, 30.5L, 2.0L, 0.3677149619770576197880557L, -0.3027246476567273059051026L},
};

template <class T> quaternion<T> eval(std::string_view fn, int n, const quaternion<T>& q) {
    if (fn == "tgamma") return quat::tgamma(q);
    if (fn == "lgamma") return quat::lgamma(q);
    if (fn == "digamma") return quat::digamma(q);
    if (fn == "erf") return quat::erf(q);
    if (fn == "erfc") return quat::erfc(q);
    if (fn == "zeta") return quat::riemann_zeta(q);
    return quat::cyl_bessel_j(n, q);
}

}  // namespace

template <class T> void reference_values() {
    for (const auto& r : refs) {
        const T a = static_cast<T>(r.re), b = static_cast<T>(r.im);
        const T wr = static_cast<T>(r.wre), wi = static_cast<T>(r.wim);
        // On the i axis, and on an arbitrary axis with the same complex plane.
        const quaternion<T> q_i(a, b), q_u = on_axis(a, b);
        const quaternion<T> want_i(wr, wi), want_u = on_axis(wr, wi);
        const bool ok_i = close(eval<T>(r.fn, r.n, q_i), want_i, rtol<T>());
        const bool ok_u = close(eval<T>(r.fn, r.n, q_u), want_u, rtol<T>());
        if (!ok_i || !ok_u) std::fprintf(stderr, "  %.*s(n=%d) at (%Lg,%Lg)\n", int(r.fn.size()), r.fn.data(), r.n, r.re, r.im);
        CHECK(ok_i);
        CHECK(ok_u);
        // Conjugate symmetry (the axis -u).
        if (r.fn != "lgamma") CHECK_CLOSE(eval<T>(r.fn, r.n, conj(q_u)), conj(want_u), rtol<T>());
    }
}
TEST(reference_values_mpmath) { FOR_EACH_TYPE(reference_values); }

// Inverse functions against std::complex along an arbitrary axis.
template <class T, class QF, class CF>
void vs_complex(QF qf, CF cf, T a, T b, T tol) {
    using W = quat::detail::sf::work_t<T>;
    const std::complex<W> z = cf(std::complex<W>(static_cast<W>(a), static_cast<W>(b)));
    const quaternion<T> want = on_axis(static_cast<T>(z.real()), static_cast<T>(z.imag()));
    const quaternion<T> got = qf(on_axis(a, b));
    if (!close(got, want, tol))
        std::printf("  vs std::complex at (%g, %g): %.3g eps, tolerance %.3g eps\n", static_cast<double>(a), static_cast<double>(b),
                    static_cast<double>(abs(got - want) / abs(want) / std::numeric_limits<T>::epsilon()),
                    static_cast<double>(tol / std::numeric_limits<T>::epsilon()));
    CHECK_CLOSE(got, want, tol);
}

template <class T> void inverse_functions() {
    using Q = quaternion<T>;
    const T tol = rtol<T>(64);
#if defined(_LIBCPP_VERSION) && defined(__APPLE__)
    // libc++'s complex inverse functions take log(abs(w)) with |w| near 1, and macOS's libm
    // hypot is not correctly rounded, so near 0 the reference itself is off by up to 387 eps
    // (asinh at 0.001 + 0.002i; the quaternion result agrees with the long double reference).
    const T ref_tol = std::is_same_v<T, float> ? tol : T(1e-12);
#elif defined(_LIBCPP_VERSION)
    // libc++ evaluates the complex<long double> inverse functions with double-precision
    // division, so as a reference it is only good to about 1e-13.
    const T ref_tol = std::is_same_v<T, long double> ? T(1e-12) : tol;
#else
    const T ref_tol = tol;
#endif
    const T pts[][2] = {{T(0.25), T(0.5)}, {T(-0.75), T(1.5)}, {T(2.5), T(-0.25)}, {T(-3), T(-2)}, {T(0.5), T(10)}, {T(1e-3), T(2e-3)}};
    for (auto& p : pts) {
#define VS(fn) vs_complex<T>([](const Q& q) { return quat::fn(q); }, [](const auto& z) { return std::fn(z); }, p[0], p[1], ref_tol)
        VS(asin); VS(acos); VS(atan); VS(asinh); VS(acosh); VS(atanh);
#undef VS
        const Q q = on_axis(p[0], p[1]);
        const T rt = rtol<T>(256);
        CHECK_CLOSE_S(sin(quat::asin(q)), q, rt, T(1));
        CHECK_CLOSE_S(cos(quat::acos(q)), q, rt, T(1));
        CHECK_CLOSE_S(tan(quat::atan(q)), q, rt, T(1));
        CHECK_CLOSE_S(sinh(quat::asinh(q)), q, rt, T(1));
        CHECK_CLOSE_S(cosh(quat::acosh(q)), q, rt, T(1));
        CHECK_CLOSE_S(tanh(quat::atanh(q)), q, rt, T(1));
        // The results stay in the plane span{1, u} of q.
        CHECK_CLOSE(quat::asin(q) * q, q * quat::asin(q), rt);
    }
    // Real axis: inside the domain the results are real, as in <cmath>.
    for (T x : {T(-0.9), T(-0.25), T(0), T(0.5), T(0.99)}) {
        CHECK(qtest::near(quat::asin(Q(x)).real(), std::asin(x), tol) && quat::asin(Q(x)).unreal() == Q());
        CHECK(qtest::near(quat::acos(Q(x)).real(), std::acos(x), tol));
        CHECK(qtest::near(quat::atanh(Q(x)).real(), std::atanh(x), tol));
        CHECK(qtest::near(quat::atan(Q(x * 7)).real(), std::atan(x * 7), tol));
        CHECK(qtest::near(quat::asinh(Q(x * 7)).real(), std::asinh(x * 7), tol));
        CHECK(qtest::near(quat::acosh(Q(T(2.5) - x)).real(), std::acosh(T(2.5) - x), tol));
    }
    // Outside the domain the cut is approached from above (imag +0), as std::complex
    // does for (x, +0), and the imaginary part goes on the i axis.
    const auto c = std::asin(std::complex<T>(T(2), T(0)));
    CHECK_CLOSE(quat::asin(Q(2)), Q(c.real(), c.imag()), tol);
    const auto h = std::acosh(std::complex<T>(T(-2), T(0)));
    CHECK_CLOSE(quat::acosh(Q(-2)), Q(h.real(), h.imag()), tol);
}
TEST(inverse_trig_and_hyperbolic) { FOR_EACH_TYPE(inverse_functions); }

template <class T> void elementary() {
    using Q = quaternion<T>;
    const T tol = rtol<T>(64);
    for (T x : {T(-3), T(-0.5), T(-1e-9), T(0), T(1e-9), T(0.75), T(5)}) {
        CHECK(qtest::near(quat::expm1(Q(x)).real(), std::expm1(x), tol) && quat::expm1(Q(x)).unreal() == Q());
        if (x > -1) CHECK(qtest::near(quat::log1p(Q(x)).real(), std::log1p(x), tol) && quat::log1p(Q(x)).unreal() == Q());
        if (x >= 0) CHECK(qtest::near(quat::cbrt(Q(x)).real(), std::cbrt(x), tol));
    }
    for (auto [a, b] : {std::pair{T(0.5), T(1.5)}, {T(-2), T(0.25)}, {T(3), T(-4)}, {T(-0.25), T(-0.5)}}) {
        const Q q = on_axis(a, b);
        CHECK_CLOSE(quat::expm1(q), exp(q) - T(1), rtol<T>(256));
        CHECK_CLOSE(quat::log1p(q), log(q + T(1)), rtol<T>(256));
        CHECK_CLOSE(pow(quat::cbrt(q), 3), q, rtol<T>(256));
        CHECK_CLOSE(quat::cbrt(q), pow(q, T(1) / 3), rtol<T>(256));
    }
    // Small arguments keep full relative accuracy (no cancellation).
    const Q tiny = on_axis(T(1e-12), T(3e-12));
    CHECK_CLOSE(quat::expm1(tiny), tiny + tiny * tiny / T(2), rtol<T>(16));
    CHECK_CLOSE(quat::log1p(tiny), tiny - tiny * tiny / T(2), rtol<T>(16));
    const T nl = std::log1p(T(-0.999999));
    CHECK(qtest::near(quat::log1p(Q(T(-0.999999))).real(), nl, tol));
    // Negative reals: principal cube root, like std::pow(complex, 1/3).
    CHECK_CLOSE(quat::cbrt(Q(-8)), Q(1, std::sqrt(T(3))), rtol<T>(64));
}
TEST(expm1_log1p_cbrt) { FOR_EACH_TYPE(elementary); }

template <class T> void gamma_family() {
    using Q = quaternion<T>;
    const T tol = rtol<T>();
    // Real axis against <cmath>.
    for (T x : {T(0.125), T(0.5), T(1), T(2.5), T(7), T(20.25), T(-0.5), T(-2.75), T(-10.5)}) {
        const Q g = quat::tgamma(Q(x));
        CHECK(qtest::near(g.real(), std::tgamma(x), tol) && g.unreal() == Q());
        CHECK(qtest::near(quat::lgamma(Q(x)).real(), std::lgamma(x), tol));
    }
    CHECK(qtest::near(quat::tgamma(Q(T(30.5))).real(), std::tgamma(T(30.5)), rtol<T>(20000)));
    CHECK(quat::tgamma(Q(5)) == Q(24) || qtest::near(quat::tgamma(Q(5)).real(), T(24), tol));
    CHECK(std::isinf(quat::tgamma(Q(0)).real()) && std::isnan(quat::tgamma(Q(-3)).real()));
    CHECK(std::isinf(quat::lgamma(Q(-2)).real()));
    // lgamma on the negative axis (limit from above): Re = log|Gamma|, Im = -pi * (number of negative factors), on i.
    CHECK_CLOSE(quat::lgamma(Q(T(-2.5))), Q(std::lgamma(T(-2.5)), -3 * quat::detail::pi_v<T>), tol);
    // Digamma: known values and recurrence.
    const T g = quat::detail::egamma_v<T>;
    CHECK(qtest::near(quat::digamma(Q(1)).real(), -g, tol));
    CHECK(qtest::near(quat::digamma(Q(T(0.5))).real(), -g - 2 * quat::detail::ln2_v<T>, tol));
    CHECK(std::isnan(quat::digamma(Q(-1)).real()));
    for (auto [a, b] : {std::pair{T(0.3), T(0.7)}, {T(-2.6), T(1.1)}, {T(4.5), T(-3)}, {T(-0.5), T(12)}, {T(25), T(5)}}) {
        const Q q = on_axis(a, b);
        const Q G = quat::tgamma(q);
        CHECK_CLOSE(quat::tgamma(q + T(1)), q * G, tol);                          // recurrence
        CHECK_CLOSE(G * quat::tgamma(T(1) - q),                                     // reflection
                    quat::detail::pi_v<T> * inverse(sin(quat::detail::pi_v<T> * q)), tol);
        CHECK_CLOSE(exp(quat::lgamma(q)), G, rtol<T>(20000));
        CHECK_CLOSE(quat::lgamma(q + T(1)), quat::lgamma(q) + log(q), tol);        // continuous branch
        CHECK_CLOSE(quat::digamma(q + T(1)), quat::digamma(q) + inverse(q), tol);
        CHECK_CLOSE(quat::digamma(T(1) - q) - quat::digamma(q),
                    quat::detail::pi_v<T> * cos(quat::detail::pi_v<T> * q) * inverse(sin(quat::detail::pi_v<T> * q)), tol);
        if (std::is_same_v<T, float> && a > 10) continue;  // Gamma(2q), Gamma(q+b) overflow float
        // Legendre duplication: Gamma(q) Gamma(q+1/2) = 2^(1-2q) sqrt(pi) Gamma(2q).
        CHECK_CLOSE(G * quat::tgamma(q + T(0.5)),
                    pow(T(2), T(1) - T(2) * q) * std::sqrt(quat::detail::pi_v<T>) * quat::tgamma(T(2) * q), rtol<T>(20000));
        // Beta.
        CHECK_CLOSE(quat::beta(q, T(1.75)), G * quat::tgamma(Q(T(1.75))) * inverse(quat::tgamma(q + T(1.75))), rtol<T>(20000));
        CHECK(quat::beta(T(1.75), q) == quat::beta(q, T(1.75)));
    }
    CHECK(qtest::near(quat::beta(Q(2), T(3)).real(), T(1) / 12, tol));
}
TEST(gamma_digamma_beta) { FOR_EACH_TYPE(gamma_family); }

template <class T> void error_functions() {
    using Q = quaternion<T>;
    const T tol = rtol<T>();
    for (T x : {T(-6), T(-2), T(-0.6), T(-0.5), T(-0.1), T(0), T(1e-8), T(0.3), T(0.5), T(1), T(3), T(8)}) {
        const Q e = quat::erf(Q(x)), c = quat::erfc(Q(x));
        CHECK(qtest::near(e.real(), std::erf(x), tol) && e.unreal() == Q());
        CHECK(qtest::near(c.real(), std::erfc(x), tol) && c.unreal() == Q());
    }
    CHECK(qtest::near(quat::erfc(Q(T(9))).real(), std::erfc(T(9)), rtol<T>(20000)));  // tiny tail keeps relative accuracy
    CHECK(quat::erf(Q(std::numeric_limits<T>::infinity())) == Q(1));
    CHECK(quat::erfc(Q(-std::numeric_limits<T>::infinity())) == Q(2));
    for (auto [a, b] : {std::pair{T(0.2), T(0.3)}, {T(-1.5), T(2)}, {T(3), T(-1)}, {T(0.45), T(6)}, {T(-4), T(-4.5)}, {T(0.55), T(0.1)}}) {
        const Q q = on_axis(a, b);
        const T sc = abs(quat::erf(q));  // erf + erfc = 1 can only hold relative to |erf|
        CHECK_CLOSE_S(quat::erf(q) + quat::erfc(q), Q(1), tol, sc);
        CHECK_CLOSE(quat::erf(-q), -quat::erf(q), tol);
        CHECK_CLOSE_S(quat::erfc(-q), T(2) - quat::erfc(q), tol, sc);
    }
    // erf(i y) = i erfi(y) = i * 2/sqrt(pi) * integral_0^y exp(t^2) dt; check with a series.
    {
        const T y = T(1.25);
        T term = y, sum = y;
        for (int n = 1; n < 200; ++n) { term *= y * y / T(n); sum += term / T(2 * n + 1); }
        CHECK_CLOSE(quat::erf(Q(0, 0, y, 0)), Q(0, 0, sum * T(2) / std::sqrt(quat::detail::pi_v<T>), 0), tol);
    }
}
TEST(erf_erfc) { FOR_EACH_TYPE(error_functions); }

template <class T> void zeta() {
    using Q = quaternion<T>;
    const T tol = rtol<T>();
    const T pi = quat::detail::pi_v<T>;
    CHECK(qtest::near(quat::riemann_zeta(Q(2)).real(), pi * pi / 6, tol));
    CHECK(qtest::near(quat::riemann_zeta(Q(4)).real(), pi * pi * pi * pi / 90, tol));
    CHECK(qtest::near(quat::riemann_zeta(Q(0)).real(), T(-0.5), tol));
    CHECK(qtest::near(quat::riemann_zeta(Q(-1)).real(), T(-1) / 12, tol));
    CHECK(qtest::near(quat::riemann_zeta(Q(-3)).real(), T(1) / 120, tol));
    CHECK(quat::riemann_zeta(Q(-2)) == Q() && quat::riemann_zeta(Q(-20)) == Q());
    CHECK(std::isinf(quat::riemann_zeta(Q(1)).real()));
    CHECK(qtest::near(quat::riemann_zeta(Q(60)).real(), T(1), tol));
    // First nontrivial zero 1/2 + 14.134725141734693790457... i.
    const Q rho = on_axis(T(0.5), static_cast<T>(14.134725141734693790457251983562L));
    CHECK(abs(quat::riemann_zeta(rho)) < rtol<T>(20000));
    // Near the pole: zeta(s) = 1/(s-1) + gamma_0 - gamma_1 (s-1) + O((s-1)^2) (Stieltjes constants).
    const Q s = on_axis(T(1) + T(1e-6), T(1e-6));
    const T gamma1 = static_cast<T>(-0.0728158454836767248605863758749547L);
    CHECK_CLOSE(quat::riemann_zeta(s), inverse(s - T(1)) + quat::detail::egamma_v<T> - gamma1 * (s - T(1)), rtol<T>(20000));
    // Functional equation, checked on an arbitrary axis.
    for (auto [a, b] : {std::pair{T(0.3), T(2)}, {T(-1.5), T(0.75)}, {T(2.5), T(-7)}}) {
        const Q q = on_axis(a, b);
        const Q rhs = pow(T(2), q) * pow(pi, q - T(1)) * sin(pi * q / T(2)) * quat::tgamma(T(1) - q) *
                      quat::riemann_zeta(T(1) - q);
        CHECK_CLOSE(quat::riemann_zeta(q), rhs, rtol<T>(20000));
    }
}
TEST(riemann_zeta) { FOR_EACH_TYPE(zeta); }

template <class T> void bessel() {
    using Q = quaternion<T>;
    const T tol = rtol<T>(20000);
    CHECK(quat::cyl_bessel_j(0, Q()) == Q(1) && quat::cyl_bessel_j(3, Q()) == Q());
    for (auto [a, b] : {std::pair{T(0.5), T(0.25)}, {T(4), T(1)}, {T(-12), T(3)}, {T(25), T(-0.5)}, {T(0.1), T(9)}}) {
        const Q q = on_axis(a, b);
        for (int n : {1, 2, 5, 9}) {
            // J_{n-1} + J_{n+1} = (2n/q) J_n
            CHECK_CLOSE(quat::cyl_bessel_j(n - 1, q) + quat::cyl_bessel_j(n + 1, q),
                        T(2 * n) * inverse(q) * quat::cyl_bessel_j(n, q), tol);
            CHECK_CLOSE(quat::cyl_bessel_j(-n, q), (n % 2 ? T(-1) : T(1)) * quat::cyl_bessel_j(n, q), tol);
            CHECK_CLOSE(quat::cyl_bessel_j(n, -q), (n % 2 ? T(-1) : T(1)) * quat::cyl_bessel_j(n, q), tol);
        }
        // Neumann addition: J0^2 + 2 sum J_n^2 = 1.
        Q s = quat::cyl_bessel_j(0, q) * quat::cyl_bessel_j(0, q);
        for (int n = 1; n < 80; ++n) s += T(2) * quat::cyl_bessel_j(n, q) * quat::cyl_bessel_j(n, q);
        CHECK(abs(s - T(1)) < tol * std::max(T(1), abs(quat::cyl_bessel_j(0, q)) * abs(quat::cyl_bessel_j(0, q))));
    }
}
TEST(bessel_j) { FOR_EACH_TYPE(bessel); }

#if HAVE_STD_SPECIAL_MATH
template <class T> void vs_std_special() {
    using Q = quaternion<T>;
    using R = quat::detail::sf::work_t<T>;  // libstdc++'s float versions are not float-accurate
    const T tol = std::max({rtol<T>(20000), T(1e-12), T(64) * std::numeric_limits<T>::epsilon()});  // libstdc++'s own error reaches ~2e-13
    for (T x : {T(-7.5), T(-2.5), T(-0.75), T(0.25), T(0.5), T(0.999), T(1.5), T(2), T(3.5), T(12)})
        CHECK(qtest::near(quat::riemann_zeta(Q(x)).real(), static_cast<T>(std::riemann_zeta(static_cast<R>(x))), tol));
    for (auto [a, b] : {std::pair{T(0.5), T(0.5)}, {T(2), T(3)}, {T(1.25), T(7.5)}, {T(10), T(0.1)}})
        CHECK(qtest::near(quat::beta(Q(a), b).real(), static_cast<T>(std::beta(static_cast<R>(a), static_cast<R>(b))), tol));
    for (int n : {0, 1, 2, 7, 20})
        for (T x : {T(0.01), T(0.5), T(1), T(3.75), T(10), T(40)}) {
            CHECK(qtest::near(quat::cyl_bessel_j(n, Q(x)).real(), static_cast<T>(std::cyl_bessel_j(R(n), static_cast<R>(x))), tol));
            // J_n(i y) = i^n I_n(y), on an arbitrary axis.
            const T in = static_cast<T>(std::cyl_bessel_i(R(n), static_cast<R>(x)));
            const Q want = (n % 4 == 0) ? Q(in) : (n % 4 == 1) ? in * axis<T>() : (n % 4 == 2) ? Q(-in) : -in * axis<T>();
            CHECK_CLOSE(quat::cyl_bessel_j(n, x * axis<T>()), want, tol);
        }
}
TEST(against_std_special_math) { FOR_EACH_TYPE(vs_std_special); }
#endif

TEST_MAIN
