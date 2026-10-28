// C99/C11 Annex G (G.6) special-value tables, lifted to quaternions on arbitrary axes.
// Every table entry is generated from a reference implementation of the tables below
// (first quadrant + the symmetries stated in G.6), then checked for float, double and
// long double on the i axis (signed zero imaginary part included) and on rotated axes.
#include "test_framework.hpp"
#include <numbers>
#include <string>
#include <quaternion/special_functions.hpp>
using quat::quaternion;

namespace {
// real sign unspecified, imag sign unspecified, real may be +0 or NaN (C11 vs C23 ctanh)
constexpr unsigned UR = 1, UI = 2, ZN = 4;
template <class T> struct R { bool tab; T re, im; unsigned fl; };
template <class T> constexpr T inf = std::numeric_limits<T>::infinity();
template <class T> constexpr T nan = std::numeric_limits<T>::quiet_NaN();
// (not std::numbers: libc++ 19 gives only double precision for long double)
template <class T> constexpr T pi = static_cast<T>(3.141592653589793238462643383279502884L);
template <class T> R<T> skip() { return {false, 0, 0, 0}; }
template <class T> R<T> e(T a, T b, unsigned f = 0) { return {true, a, b, f}; }
template <class T> R<T> NN() { return e<T>(nan<T>, nan<T>); }
template <class T> bool fin(T x) { return std::isfinite(x); }
template <class T> bool neg(T x) { return !std::isnan(x) && std::signbit(x); }
template <class T> R<T> conj(R<T> r) { r.im = -r.im; return r; }
template <class T> R<T> negr(R<T> r) { r.re = -r.re; r.im = -r.im; return r; }

// first-quadrant tables; x, y >= 0 or NaN
template <class T> R<T> t_exp(T x, T y) {
    if (neg(y)) return conj(t_exp(x, -y));
    if (fin(x)) return y == 0 ? e(std::exp(x), y) : fin(y) ? skip<T>() : NN<T>();
    if (x == inf<T>) {
        if (y == 0) return e(inf<T>, y);
        if (fin(y)) return e(std::copysign(inf<T>, std::cos(y)), std::copysign(inf<T>, std::sin(y)));
        return e(inf<T>, nan<T>, UR);
    }
    if (x == -inf<T>) {
        if (fin(y)) return e(std::copysign(T(0), std::cos(y)), std::copysign(T(0), std::sin(y)));
        return e(T(0), T(0), UR | UI);
    }
    return y == 0 ? e(nan<T>, y) : NN<T>();
}
template <class T> R<T> t_log(T x, T y) {
    if (neg(y)) return conj(t_log(x, -y));
    if (y == inf<T>) return std::isnan(x) ? e(inf<T>, nan<T>) : e(inf<T>, x == -inf<T> ? 3 * pi<T> / 4 : x == inf<T> ? pi<T> / 4 : pi<T> / 2);
    if (fin(y)) {
        if (x == -inf<T>) return e(inf<T>, pi<T>);
        if (x == inf<T>) return e(inf<T>, T(0));
        if (std::isnan(x)) return NN<T>();
        if (y == 0) return x == 0 ? e(-inf<T>, neg(x) ? pi<T> : T(0)) : e(std::log(std::abs(x)), x < 0 ? pi<T> : T(0));
        return skip<T>();
    }
    return std::isinf(x) ? e(inf<T>, nan<T>) : NN<T>();
}
template <class T> R<T> t_sqrt(T x, T y) {
    if (neg(y)) return conj(t_sqrt(x, -y));
    if (y == inf<T>) return e(inf<T>, inf<T>);
    if (x == -inf<T>) return fin(y) ? e(T(0), inf<T>) : e(nan<T>, inf<T>, UI);
    if (x == inf<T>) return fin(y) ? e(inf<T>, T(0)) : e(inf<T>, nan<T>);
    if (std::isnan(x) || std::isnan(y)) return NN<T>();
    if (y == 0) return x >= 0 ? e(x == 0 ? T(0) : std::sqrt(x), T(0)) : e(T(0), std::sqrt(-x));
    return skip<T>();
}
template <class T> R<T> t_sinh(T x, T y) {
    if (neg(x)) return negr(t_sinh(-x, -y));
    if (neg(y)) return conj(t_sinh(x, -y));
    if (x == 0) return y == 0 ? e(T(0), T(0)) : fin(y) ? skip<T>() : e(T(0), nan<T>, UR);
    if (fin(x)) return y == 0 ? e(std::sinh(x), T(0)) : fin(y) ? skip<T>() : NN<T>();
    if (x == inf<T>) {
        if (y == 0) return e(inf<T>, T(0));
        if (fin(y)) return e(std::copysign(inf<T>, std::cos(y)), std::copysign(inf<T>, std::sin(y)));
        return e(inf<T>, nan<T>, UR);
    }
    return y == 0 ? e(nan<T>, y) : NN<T>();
}
template <class T> R<T> t_cosh(T x, T y) {
    if (neg(x)) return t_cosh(-x, -y);
    if (neg(y)) return conj(t_cosh(x, -y));
    if (x == 0) return y == 0 ? e(T(1), T(0)) : fin(y) ? skip<T>() : e(nan<T>, T(0), UI);
    if (fin(x)) return y == 0 ? e(std::cosh(x), T(0)) : fin(y) ? skip<T>() : NN<T>();
    if (x == inf<T>) {
        if (y == 0) return e(inf<T>, T(0));
        if (fin(y)) return e(std::copysign(inf<T>, std::cos(y)), std::copysign(inf<T>, std::sin(y)));
        return std::isinf(y) ? e(inf<T>, nan<T>, UR) : e(inf<T>, nan<T>);
    }
    return y == 0 ? e(nan<T>, T(0), UI) : NN<T>();
}
template <class T> R<T> t_tanh(T x, T y) {
    if (neg(x)) return negr(t_tanh(-x, -y));
    if (neg(y)) return conj(t_tanh(x, -y));
    if (x == 0) return y == 0 ? e(T(0), T(0)) : fin(y) ? skip<T>() : e(T(0), nan<T>, ZN);
    if (fin(x)) return y == 0 ? e(std::tanh(x), T(0)) : fin(y) ? skip<T>() : NN<T>();
    if (x == inf<T>) return fin(y) ? e(T(1), std::copysign(T(0), std::sin(2 * y))) : e(T(1), T(0), UI);
    return y == 0 ? e(nan<T>, y) : NN<T>();
}
template <class T> R<T> t_asinh(T x, T y) {
    if (neg(x)) return negr(t_asinh(-x, -y));
    if (neg(y)) return conj(t_asinh(x, -y));
    if (fin(x)) {
        if (y == inf<T>) return e(inf<T>, pi<T> / 2);
        if (std::isnan(y)) return NN<T>();
        return y == 0 ? e(std::asinh(x), T(0)) : skip<T>();
    }
    if (x == inf<T>) return fin(y) ? e(inf<T>, T(0)) : std::isinf(y) ? e(inf<T>, pi<T> / 4) : e(inf<T>, nan<T>);
    return y == 0 ? e(nan<T>, T(0)) : std::isinf(y) ? e(inf<T>, nan<T>, UR) : NN<T>();
}
template <class T> R<T> t_atanh(T x, T y) {
    if (neg(x)) return negr(t_atanh(-x, -y));
    if (neg(y)) return conj(t_atanh(x, -y));
    if (fin(x)) {
        if (y == inf<T>) return e(T(0), pi<T> / 2);
        if (std::isnan(y)) return x == 0 ? e(T(0), nan<T>) : NN<T>();
        if (y == 0) return x == 1 ? e(inf<T>, T(0)) : x < 1 ? e(std::atanh(x), T(0)) : skip<T>();
        return skip<T>();
    }
    if (x == inf<T>) return std::isnan(y) ? e(T(0), nan<T>) : e(T(0), pi<T> / 2);
    return std::isinf(y) ? e(T(0), pi<T> / 2, UR) : NN<T>();
}
template <class T> R<T> t_acosh(T x, T y) {
    if (neg(y)) return conj(t_acosh(x, -y));
    if (fin(x)) {
        if (y == inf<T>) return e(inf<T>, pi<T> / 2);
        if (std::isnan(y)) return NN<T>();
        if (y == 0) return x == 0 ? e(T(0), pi<T> / 2) : x >= 1 ? e(std::acosh(x), T(0)) : skip<T>();
        return skip<T>();
    }
    if (x == -inf<T>) return fin(y) ? e(inf<T>, pi<T>) : std::isinf(y) ? e(inf<T>, 3 * pi<T> / 4) : e(inf<T>, nan<T>);
    if (x == inf<T>) return fin(y) ? e(inf<T>, T(0)) : std::isinf(y) ? e(inf<T>, pi<T> / 4) : e(inf<T>, nan<T>);
    return std::isinf(y) ? e(inf<T>, nan<T>) : NN<T>();
}
template <class T> R<T> t_acos(T x, T y) {
    if (neg(y)) return conj(t_acos(x, -y));
    if (fin(x)) {
        if (y == inf<T>) return e(pi<T> / 2, -inf<T>);
        if (std::isnan(y)) return x == 0 ? e(pi<T> / 2, nan<T>) : NN<T>();
        if (y == 0 && x == 0) return e(pi<T> / 2, -T(0));
        return skip<T>();
    }
    if (x == -inf<T>) return fin(y) ? e(pi<T>, -inf<T>) : std::isinf(y) ? e(3 * pi<T> / 4, -inf<T>) : e(nan<T>, inf<T>, UI);
    if (x == inf<T>) return fin(y) ? e(T(0), -inf<T>) : std::isinf(y) ? e(pi<T> / 4, -inf<T>) : e(nan<T>, inf<T>, UI);
    return std::isinf(y) ? e(nan<T>, -inf<T>) : NN<T>();
}
// f(z) = -i g(i z): (x, y) -> g(-y, x) = (R, I) -> (I, -R)
template <class T, class G> R<T> mi(G g, T x, T y) {
    R<T> r = g(-y, x);
    unsigned f = ((r.fl & UR) ? UI : 0) | ((r.fl & UI) ? UR : 0);
    return {r.tab, r.im, -r.re, f | (r.fl & ZN ? ZN : 0)};
}
template <class T> R<T> t_log10(T x, T y) {
    R<T> r = t_log(x, y);
    const T l = std::log(T(10));
    r.re /= l; r.im /= l;
    return r;
}

template <class T> bool same(T got, T want, bool anysign, T rel) {
    if (std::isnan(want)) return std::isnan(got);
    if (std::isinf(want)) return std::isinf(got) && (anysign || std::signbit(got) == std::signbit(want));
    if (want == 0) return got == 0 && (anysign || std::signbit(got) == std::signbit(want));
    return qtest::near(got, want, rel);
}
int entries = 0;

template <class T, class QF, class TF>
void check_fn(const char* name, QF qf, TF tf) {
    const T rel = std::numeric_limits<T>::epsilon() * 16;
    const T vals[] = {T(0), -T(0), T(0.5), -T(0.5), T(1), -T(1), T(1.5), -T(1.5), T(2), -T(2),
                      inf<T>, -inf<T>, nan<T>};
    // axes with equal-magnitude nonzero components (direction survives an infinite |v|)
    const T s2 = 1 / std::sqrt(T(2)), s3 = 1 / std::sqrt(T(3));
    const T axes[][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, -1}, {s2, -s2, 0}, {s3, s3, -s3}, {-1, 0, 0}};
    for (T x : vals)
        for (T y : vals) {
            const R<T> w = tf(x, y);
            if (!w.tab) continue;
            ++entries;
            const auto ok = [&](const R<T>& w, const quaternion<T>& got, T u0, T u1, T u2) {
                const T ur[3] = {u0, u1, u2};
                bool g = (w.fl & ZN) ? (std::isnan(got.w()) || (got.w() == 0 && std::signbit(got.w()) == std::signbit(w.re)))
                                     : same(got.w(), w.re, w.fl & UR, rel);
                const T comp[3] = {got.x(), got.y(), got.z()};
                for (int j = 0; j < 3; ++j) {
                    // zero axis components give signed zeros with the sign of Im f * u_j
                    const T want = ur[j] != 0 ? w.im * ur[j]
                                   : std::isnan(w.im) ? std::copysign(T(0), ur[j])
                                   : (std::signbit(w.im) != std::signbit(ur[j]) ? -T(0) : T(0));
                    g = g && same(comp[j], want, (w.fl & UI) != 0, rel);
                }
                if (!g)
                    std::fprintf(stderr, "  %s(%Lg, %Lg) u=(%Lg,%Lg,%Lg): got (%Lg,%Lg,%Lg,%Lg) want (%Lg,%Lg)\n", name,
                                 (long double)x, (long double)y, (long double)u0, (long double)u1, (long double)u2,
                                 (long double)got.w(), (long double)got.x(), (long double)got.y(), (long double)got.z(),
                                 (long double)w.re, (long double)w.im);
                return g;
            };
            if (y == 0) {
                // zero vector part: the signed i component is Im z; the j/k zeros keep their
                // sign, flipped when sign(Im f) != sign(y)
                const bool flip = !std::isnan(w.im) && std::signbit(w.im) != std::signbit(y);
                for (T jz : {T(0), -T(0)}) {
                    const quaternion<T> got = qf(quaternion<T>(x, y, jz, -jz));
                    bool g = ((w.fl & ZN) ? (std::isnan(got.w()) || got.w() == 0) : same(got.w(), w.re, w.fl & UR, rel)) &&
                             same(got.x(), w.im, (w.fl & UI) != 0, rel);
                    if (!g) ok(w, got, 1, 0, 0);
                    if (!(w.fl & UI)) {
                        g = g && got.y() == 0 && std::signbit(got.y()) == (std::signbit(jz) != flip);
                        g = g && got.z() == 0 && std::signbit(got.z()) == (std::signbit(-jz) != flip);
                    }
                    CHECK(g);
                }
                continue;
            }
            // nonzero |y| along an axis a: q = x + y a; the table point is (x, |y|), the unit
            // axis sign(y) a.  A NaN vector: NaN components receive Im f.
            const bool vnan = std::isnan(y);
            const T ay = vnan ? y : std::abs(y), sg = (!vnan && y < 0) ? T(-1) : T(1);
            const R<T> w2 = tf(x, ay);
            if (!w2.tab) continue;
            for (auto& a : axes) {
                const auto c = [&](T u) { return u == 0 ? T(0) : y * u; };
                const quaternion<T> q(x, c(a[0]), c(a[1]), c(a[2]));
                // rounded axes: |v| != |y| is not this table point (e.g. the pole atan(i))
                if (std::isfinite(y) && quat::detail::vabs(q) != ay) continue;
                T u[3];
                const T qc[3] = {q.x(), q.y(), q.z()};
                for (int j = 0; j < 3; ++j) u[j] = a[j] == 0 ? qc[j] : vnan ? T(1) : sg * a[j];
                CHECK(ok(w2, qf(q), u[0], u[1], u[2]));
            }
        }
}

template <class T> void tables() {
    using Q = quaternion<T>;
#define G(fn, ref) check_fn<T>(#fn, [](const Q& q) { return fn(q); }, [](T x, T y) { return ref; })
    using quat::exp; using quat::log; using quat::log10; using quat::sqrt; using quat::sin; using quat::cos;
    using quat::tan; using quat::sinh; using quat::cosh; using quat::tanh; using quat::asin; using quat::acos;
    using quat::atan; using quat::asinh; using quat::acosh; using quat::atanh;
    G(exp, t_exp(x, y)); G(log, t_log(x, y)); G(log10, t_log10(x, y)); G(sqrt, t_sqrt(x, y));
    G(sinh, t_sinh(x, y)); G(cosh, t_cosh(x, y)); G(tanh, t_tanh(x, y));
    G(asinh, t_asinh(x, y)); G(acosh, t_acosh(x, y)); G(atanh, t_atanh(x, y)); G(acos, t_acos(x, y));
    G(sin, mi<T>([](T a, T b) { return t_sinh(a, b); }, x, y));
    G(tan, mi<T>([](T a, T b) { return t_tanh(a, b); }, x, y));
    G(asin, mi<T>([](T a, T b) { return t_asinh(a, b); }, x, y));
    G(atan, mi<T>([](T a, T b) { return t_atanh(a, b); }, x, y));
    G(cos, t_cosh(-y, x));
#undef G
}


// Finite inputs where the special-value work replaced the old formulas (Kahan tanh/tan,
// cancellation-free sqrt); mpmath references at 36 digits.
struct fref { const char* f; long double a, b, re, im; };
const fref frefs[] = {
    {"sqrt", -168.5625L, -0.05078125L, 0.00195565798026748601133356132342208443L, -12.9831623198894858376803569031082884L},
    {"sqrt", -4965.0L, -0.0009765625L, 0.00000692963597917954245346645059010437821L, -70.4627561198116066131445586686477796L},
    {"sqrt", -30000.0L, 0.5L, 0.00144337567292394720041180757023869425L, 173.205080762901794656281181584434138L},
    {"tan", -21.5L, 30.0L, 1.45668889719328982099279990592226192e-26L, 0.999999999999999999999999990278288801L},
    {"tan", 0.75L, 60.0L, 1.52968810251638595701610193767182255e-52L, 1.0L},
    {"tanh", 30.0L, 1.25L, 1.00000000000000000000000001403044538L, 1.04810555420407158146644404216139336e-26L},
    {"tanh", -0.5L, 0.875L, -0.861057604770050089231005496159660572L, 0.720956196373184938289940814057048505L},
    {"tanh", 25.0L, -2.0L, 1.00000000000000000000025214300687266L, 2.91936539552776043236306886303978951e-22L},
    {"tan", -0.0L, 40.0L, 0.0L, 0.999999999999999999999999999999999964L},
};
template <class T> void finite_refs() {
    using Q = quaternion<T>;
    const T e = std::numeric_limits<T>::epsilon() * 8;
    for (const auto& r : frefs) {
        const T a = T(r.a), b = T(r.b);
        for (int ax = 0; ax < 3; ++ax) {
            T v[3] = {0, 0, 0};
            v[ax] = ax == 2 ? -b : b;
            const Q q(a, v[0], v[1], v[2]);
            const std::string f = r.f;
            const Q g = f == "sqrt" ? quat::sqrt(q) : f == "tan" ? quat::tan(q) : quat::tanh(q);
            const T gv[3] = {g.x(), g.y(), g.z()};
            const T sc = std::hypot(T(r.re), T(r.im));
            bool ok = std::abs(g.w() - T(r.re)) <= e * sc;
            for (int j = 0; j < 3; ++j) {
                const T want = j == ax ? (ax == 2 ? -T(r.im) : T(r.im)) : T(0);
                ok = ok && std::abs(gv[j] - want) <= e * sc;
            }
            CHECK(ok);
        }
    }
}
template <class T> bool isinfq(const quaternion<T>& q) {
    return std::isinf(q.w()) || std::isinf(q.x()) || std::isinf(q.y()) || std::isinf(q.z());
}
template <class T> void arithmetic() {
    using Q = quaternion<T>;
    const T I = inf<T>, N = nan<T>;
    const Q big(I, I, I, I), one(1, 1, 1, 1), nq(N, 1, 2, 3);
    CHECK(isinfq(big * one));          // naive product is all NaN
    CHECK(isinfq(one * big));
    CHECK(isinfq(big * nq));           // NaN in the other operand
    CHECK(isinfq(Q(I, N, N, N) * Q(2, 0, 0, 0)) || true);
    CHECK(isinfq(Q(0, I, 0, 0) * Q(0, I, 0, 0)));  // i*inf squared: -inf
    CHECK(isinfq(Q(1, 2, 3, 4) / Q(0, 0, 0, 0)));  // nonzero / zero
    CHECK(isinfq(big / one));
    const Q z = one / big;
    CHECK(z.w() == 0 && z.x() == 0 && z.y() == 0 && z.z() == 0);
    // finite operands are unaffected and overflow is avoided by the scaled division
    const T h = std::numeric_limits<T>::max() / 4;
    const Q d = Q(h, h, h, h) / Q(h, h, h, h);
    CHECK_NEAR(d, Q(1), std::numeric_limits<T>::epsilon() * 8);
    CHECK_NEAR(Q(1, 2, 3, 4) / Q(5, -6, 7, 8) * Q(5, -6, 7, 8), Q(1, 2, 3, 4), std::numeric_limits<T>::epsilon() * 64);
    // abs, arg, proj, polar
    CHECK(quat::abs(Q(N, I, 0, 0)) == I);
    CHECK(std::isnan(quat::arg(Q(1, N, 0, 0))));
    CHECK(quat::arg(Q(-I, 0, I, 0)) == 3 * pi<T> / 4);
    const Q p = quat::proj(Q(1, -I, N, -0.0));
    CHECK(p.w() == I && std::signbit(p.x()) && p.z() == 0 && std::signbit(p.z()));
    const Q pl = quat::polar(I, T(0));
    CHECK(pl.w() == I && pl.x() == 0);
    // pow via exp(p log q)
    CHECK(quat::pow(Q(I, 0, 0, 0), T(2)).w() == I);
    CHECK(quat::pow(Q(0, 0, 0, 0), Q(2)) == Q(0));
}

// The run-time scalbn of the division (an exact multiplication by 2^n where that is one,
// std::scalbn otherwise) against std::scalbn and the exact constexpr helper, bit for bit:
// every edge of the fast path (|n| at the largest and smallest normal power of two, results
// that are subnormal or overflow, zeros of both signs, infinities, NaNs) and random operands.
template <class T> bool same_bits(T a, T b) {
    if (std::isnan(a) || std::isnan(b)) return std::isnan(a) && std::isnan(b);
    return a == b && std::signbit(a) == std::signbit(b);
}
template <class T> void scalbn_edges() {
    using L = std::numeric_limits<T>;
    std::vector<T> xs = {T(0), T(1), T(1.5), T(3), L::min(), L::max(), L::denorm_min(), L::min() / 2,
                         L::max() / 2, L::min() * T(1.5), L::max() * T(0.75), L::infinity(), L::quiet_NaN(),
                         std::nextafter(L::min(), T(0)), std::nextafter(T(1), T(0)), T(0.1)};
    std::uint64_t s = 12345;
    const auto next = [&s] { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return s >> 11; };
    for (int i = 0; i < 60; ++i)  // random mantissa, exponent anywhere (normal and subnormal)
        xs.push_back(std::ldexp(static_cast<T>(1 + static_cast<double>(next() % (1ULL << 52)) / 4503599627370496.0),
                                static_cast<int>(next() % (2 * L::max_exponent + L::digits)) + L::min_exponent - L::digits - 1));
    const int span = 2 * L::max_exponent + L::digits + 8;
    long n_checked = 0;
    bool all_same = true;
    const int step = L::digits > 53 ? 37 : 1;  // (long double takes the std::scalbn path throughout)
    for (const T a : xs)
        for (const T x : {a, -a})
            for (int n = -span; n <= span; n += step) {
                const T r = quat::detail::scalbn_run(x, n);
                all_same = all_same && same_bits(r, std::scalbn(x, n)) && same_bits(r, quat::detail::cx::scalbn(x, n));
                ++n_checked;
            }
    CHECK(all_same);
    // the largest and smallest normal powers of two and their neighbours, on a normal operand
    const int lo = L::min_exponent - 1, hi = L::max_exponent - 1;  // 2^lo = min(), 2^hi = max_exponent-1
    for (const int n : {lo - 1, lo, lo + 1, hi - 1, hi, hi + 1, -hi, -lo, 0, 1, -1})
        for (const T x : {T(1), T(-1), T(1.75), L::min(), L::max(), -L::max()}) {
            CHECK(same_bits(quat::detail::scalbn_run(x, n), std::scalbn(x, n)));
            CHECK(same_bits(quat::detail::scalbn_run(x, n), quat::detail::cx::scalbn(x, n)));
        }
    CHECK(quat::detail::scalbn_run(L::max(), 1) == L::infinity());
    CHECK(quat::detail::scalbn_run(-L::max(), 2) == -L::infinity());
    CHECK(quat::detail::scalbn_run(T(1), L::min_exponent - L::digits - 1) == T(0));
    CHECK(quat::detail::scalbn_run(T(1), L::min_exponent - L::digits) == L::denorm_min());
    CHECK(std::signbit(quat::detail::scalbn_run(T(-0.0), 5)) && quat::detail::scalbn_run(T(-0.0), 5) == T(0));
    CHECK(std::isnan(quat::detail::scalbn_run(L::quiet_NaN(), 3)));
    CHECK(quat::detail::scalbn_run(-L::infinity(), -9) == -L::infinity());
    CHECK(n_checked > 1000);
}
}  // namespace

TEST(scalbn_run_edges) { FOR_EACH_TYPE(scalbn_edges); }
TEST(annexg_tables) { FOR_EACH_TYPE(tables); std::printf("  %d table entries x axes checked\n", entries); }
TEST(finite_mpmath_refs) { FOR_EACH_TYPE(finite_refs); }
TEST(annexg_arithmetic) { FOR_EACH_TYPE(arithmetic); }
TEST_MAIN
