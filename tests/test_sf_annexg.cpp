// Special values of the special functions (README, "Special functions: infinities, NaNs and
// signed zeros").  The tables restate the policy point by point: on every ray of a non-finite
// argument each component is the limit of that component of f, NaN where it has none; a NaN
// gives NaN; poles and zeros at the listed finite points of the real axis.  Each entry is
// checked for float, double and long double on the i axis (signed zero included) and on
// rotated axes, as in test_annexg.  Polynomials are checked against their own recurrence
// evaluated far out on the ray (an independent oracle).
#include "test_framework.hpp"
#include <complex>
#include <string>
#include <vector>
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


template <class T> R<T> zz() { return e<T>(0, 0, UR | UI); }
template <class T> T ci(T c) { return c == 0 ? c : std::copysign(inf<T>, c); }
template <class T> R<T> cinf(T c, T s, unsigned f = 0) { return e<T>(ci(c), ci(s), f | (c == 0 ? UR : 0) | (s == 0 ? UI : 0)); }
template <class T> T cospi_(T p) { T r = std::fmod(std::abs(p), T(2)); if (r == T(0.5) || r == T(1.5)) return 0; return std::cos(pi<T> * p); }
template <class T> T sinpi_(T p) { if (p == std::floor(p)) return 0; return std::sin(pi<T> * p); }
// cos/sin of (p pi - x)
template <class T> R<T> cphase(T p, T x) {
    const T c = x == 0 ? cospi_(p) : cospi_(p) * std::cos(x) + sinpi_(p) * std::sin(x);
    const T s = x == 0 ? sinpi_(p) : sinpi_(p) * std::cos(x) - cospi_(p) * std::sin(x);
    return cinf(c, s);
}
template <class T> bool isint(T x) { return fin(x) && x == std::floor(x); }
enum Ray { RP, RM, HP, HM, VV, DP, DM, FIN };
template <class T> Ray kind(T x, T y) {
    if (fin(x) && fin(y)) return FIN;
    if (std::isinf(y)) return fin(x) ? VV : x > 0 ? DP : DM;
    if (y == 0) return x > 0 ? RP : RM;
    return x > 0 ? HP : HM;
}
// conjugate-symmetric table: NaN rule, y <= -0 by conjugation, f(x, y, kind) for y >= +0
template <class T, class F> R<T> tab(T x, T y, bool realnan, F f) {
    if (std::isnan(x) || std::isnan(y)) return (y == 0 && realnan) ? e(nan<T>, y) : NN<T>();
    if (neg(y)) return conj(tab(x, -y, realnan, f));
    return f(x, y, kind(x, y));
}
#define TB(REALNAN, ...) tab<T>(x, y, REALNAN, [&](T x, T y, Ray k) -> R<T> { (void)x; (void)y; (void)k; __VA_ARGS__ })

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



// ---- polynomial oracle: recurrence in std::complex<long double> at t = 2^60 on the ray ------
using CL = std::complex<long double>;
template <class T> R<T> oracle(T x, T y, bool realnan, CL (*p)(unsigned, unsigned, CL), unsigned n, unsigned m) {
    if (std::isnan(x) || std::isnan(y)) return (y == 0 && realnan) ? e(nan<T>, y) : NN<T>();
    if (fin(x) && fin(y)) return skip<T>();
    if (neg(y)) return conj(oracle(x, -y, realnan, p, n, m));
    const long double t = 1152921504606846976.0L;  // 2^60
    const long double zx = std::isinf(x) ? (x > 0 ? t : -t) : static_cast<long double>(x);
    const long double zy = std::isinf(y) ? t : static_cast<long double>(y);
    const CL v = p(n, m, CL(zx, zy));
    const auto cl = [](long double c, bool& unspec) -> T {
        if (std::abs(c) > 1e9L) return std::copysign(inf<T>, static_cast<T>(c));
        if (std::abs(c) < 1e-9L) { unspec = true; return 0; }
        return static_cast<T>(c);
    };
    bool ur = false, ui = false;
    const T re = cl(v.real(), ur);
    T im = cl(v.imag(), ui);
    if (y == 0 && realnan) { im = y; ui = false; }  // real polynomial on the real axis
    return e(re, im, (ur ? UR : 0) | (ui ? UI : 0));
}
CL p_legendre(unsigned n, unsigned, CL x) {
    CL p0(1), p1 = x;
    if (n == 0) return p0;
    for (unsigned k = 1; k < n; ++k) { CL p2 = (2.0L * k + 1) * x * p1 / (k + 1.0L) - (long double)k * p0 / (k + 1.0L); p0 = p1; p1 = p2; }
    return p1;
}
CL p_hermite(unsigned n, unsigned, CL x) {
    CL p0(1), p1 = 2.0L * x;
    if (n == 0) return p0;
    for (unsigned k = 1; k < n; ++k) { CL p2 = 2.0L * x * p1 - 2.0L * k * p0; p0 = p1; p1 = p2; }
    return p1;
}
CL p_laguerre(unsigned n, unsigned m, CL x) {
    CL p0(1), p1 = 1.0L + (long double)m - x;
    if (n == 0) return p0;
    for (unsigned k = 1; k < n; ++k) { CL p2 = ((2.0L * k + 1 + m - x) * p1 - (k + (long double)m) * p0) / (k + 1.0L); p0 = p1; p1 = p2; }
    return p1;
}
CL p_cheb_t(unsigned n, unsigned, CL x) {
    CL p0(1), p1 = x;
    if (n == 0) return p0;
    for (unsigned k = 1; k < n; ++k) { CL p2 = 2.0L * x * p1 - p0; p0 = p1; p1 = p2; }
    return p1;
}
CL p_cheb_u(unsigned n, unsigned, CL x) {
    CL p0(1), p1 = 2.0L * x;
    if (n == 0) return p0;
    for (unsigned k = 1; k < n; ++k) { CL p2 = 2.0L * x * p1 - p0; p0 = p1; p1 = p2; }
    return p1;
}
// P_l^m with s = sqrt(1 - z^2): principal off the axis, +i sqrt(x^2 - 1) on it (library convention)
CL p_assoc_legendre(unsigned l, unsigned m, CL x) {
    CL s = x.imag() == 0 && std::abs(x.real()) > 1 ? CL(0, std::sqrt(x.real() * x.real() - 1)) : std::sqrt(1.0L - x * x);
    CL pmm(1);
    for (unsigned i = 1; i <= m; ++i) pmm *= s * (2.0L * i - 1);
    if (l == m) return pmm;
    CL pm1 = (2.0L * m + 1) * x * pmm;
    for (unsigned k = m + 2; k <= l; ++k) { CL p2 = ((2.0L * k - 1) * x * pm1 - (k + m - 1.0L) * pmm) / (long double)(k - m); pmm = pm1; pm1 = p2; }
    return pm1;
}

template <class T> void tables() {
    using Q = quaternion<T>;
    const T h = pi<T> / 2;
#define G(NAME, CALL, TABLE) check_fn<T>(NAME, [&](const Q& q) { return CALL; }, [&](T x, T y) -> R<T> { return TABLE; })
    // ---- gamma family
    G("tgamma", quat::tgamma(q), TB(true,
        if (k == FIN) return y == 0 && x == 0 ? e(std::copysign(inf<T>, x), y) : y == 0 && x < 0 && isint(x) ? e(nan<T>, y) : skip<T>();
        if (k == RP) return e(inf<T>, y);
        if (k == RM) return e(nan<T>, y);
        if (k == HP || k == DP) return NN<T>();
        return zz<T>();));
    G("lgamma", quat::lgamma(q), TB(false,
        if (k == FIN) return skip<T>();
        if (k == RP) return e(inf<T>, y);
        if (k == RM) return e(inf<T>, -inf<T>);
        if (k == HP || k == DP) return e(inf<T>, inf<T>);
        if (k == HM) return e(-inf<T>, -inf<T>);
        return e(-inf<T>, inf<T>);));
    const auto digamma_t = [&](T x, T y) { return TB(true,
        if (k == FIN) return y == 0 && x == 0 ? e(-std::copysign(inf<T>, x), y) : y == 0 && x < 0 && isint(x) ? e(nan<T>, y) : skip<T>();
        if (k == RP) return e(inf<T>, y);
        if (k == RM) return e(nan<T>, y);
        if (k == HP) return e(inf<T>, T(0), UI);
        if (k == HM) return e(inf<T>, nan<T>);
        if (k == VV) return e(inf<T>, h);
        if (k == DP) return e(inf<T>, h / 2);
        return e(inf<T>, 3 * h / 2);); };
    G("digamma", quat::digamma(q), digamma_t(x, y));
    G("polygamma0", quat::polygamma(0u, q), digamma_t(x, y));
    for (unsigned m : {1u, 2u, 3u})
        G(("polygamma" + std::to_string(m)).c_str(), quat::polygamma(m, q), TB(true,
            if (k == FIN) {
                if (y != 0) return skip<T>();
                if (x == 0) return e(m % 2 ? inf<T> : -std::copysign(inf<T>, x), y);
                if (x < 0 && isint(x)) return e(m % 2 ? inf<T> : nan<T>, y);
                return skip<T>();
            }
            if (k == RP) return e(T(0), y, UR);
            if (k == RM) return e(nan<T>, y);
            if (k == HM) return NN<T>();
            return zz<T>();));
    G("trigamma", quat::trigamma(q), TB(true,
        if (k == FIN) return y == 0 && (x == 0 || (x < 0 && isint(x))) ? e(inf<T>, y) : skip<T>();
        if (k == RP) return e(T(0), y, UR);
        if (k == RM) return e(nan<T>, y);
        if (k == HM) return NN<T>();
        return zz<T>();));
    // ---- error functions, zeta
    G("erf", quat::erf(q), TB(true,
        if (k == FIN) return y == 0 && x == 0 ? e(x, y) : skip<T>();
        if (k == RP) return e(T(1), y);
        if (k == RM) return e(T(-1), y);
        if (k == HP || k == DP) return e(T(1), T(0), UI);
        if (k == HM || k == DM) return e(T(-1), T(0), UI);
        return x == 0 ? e(x, inf<T>) : NN<T>();));
    G("erfc", quat::erfc(q), TB(true,
        if (k == FIN) return skip<T>();
        if (k == RP) return e(T(0), y, UR);
        if (k == RM) return e(T(2), y);
        if (k == HP || k == DP) return zz<T>();
        if (k == HM || k == DM) return e(T(2), T(0), UI);
        return x == 0 ? e(T(1), -inf<T>) : NN<T>();));
    G("riemann_zeta", quat::riemann_zeta(q), TB(true,
        if (k == FIN) return skip<T>();
        if (k == RP) return e(T(1), y);
        if (k == RM) return e(nan<T>, y);
        if (k == HP || k == DP) return e(T(1), T(0), UI);
        return NN<T>();));
    for (T a : {T(0.5), T(1), T(2), T(-0.5)})
        G(("hurwitz_zeta a=" + std::to_string((double)a)).c_str(), quat::hurwitz_zeta(q, a), TB(a > 0,
            if (k == FIN) return skip<T>();
            if (!(a > 0)) return NN<T>();
            if (k == RM) return e(nan<T>, y);
            if (k != RP && k != HP && k != DP) return NN<T>();
            if (a == 1) return k == RP ? e(T(1), y) : e(T(1), T(0), UI);
            if (a > 1) return k == RP ? e(T(0), y, UR) : zz<T>();
            if (k == RP) return e(inf<T>, y);
            if (k == DP) return NN<T>();
            return cinf(std::cos(y * std::log(a)), -std::sin(y * std::log(a)));));
    // ---- Bessel functions
    for (int n : {0, 1, 2, -3})
        G(("cyl_bessel_j n=" + std::to_string(n)).c_str(), quat::cyl_bessel_j(n, q), TB(true,
            if (k == FIN) {
                if (y != 0 || x != 0 || n == 0) return skip<T>();
                const bool negz = std::signbit(x) && (n % 2 != 0);
                return e((negz != (n < 0 && n % 2 != 0)) ? -T(0) : T(0), y);
            }
            if (k == RP || k == RM) return e(T(0), y, UR);
            if (k == HP || k == HM) return zz<T>();
            if (k == VV) return cphase(T(n) / 2, x);
            return NN<T>();));
    for (T nu : {T(0.5), T(2.25), T(1)})
        G(("cyl_bessel_j nu=" + std::to_string((double)nu)).c_str(), quat::cyl_bessel_j(nu, q), TB(nu == 1,
            if (k == FIN) return nu == 1 && y == 0 && x == 0 ? e(x, y) : skip<T>();
            if (k == RP || k == RM) return e(T(0), y, UR | UI);
            if (k == HP || k == HM) return zz<T>();
            if (k == VV) return cphase(nu / 2, x);
            return NN<T>();));
    for (T nu : {T(0), T(1.5)}) {
        G(("cyl_neumann nu=" + std::to_string((double)nu)).c_str(), quat::cyl_neumann(nu, q), TB(false,
            if (k == FIN) return skip<T>();
            if (k == RP || k == RM) return e(T(0), y, UR | UI);
            if (k == HP || k == HM) return zz<T>();
            if (k == VV) return cphase((nu + 1) / 2, x);
            return NN<T>();));
        G(("cyl_bessel_k nu=" + std::to_string((double)nu)).c_str(), quat::cyl_bessel_k(nu, q), TB(false,
            if (k == FIN) return skip<T>();
            if (k == RP) return e(T(0), y, UR);
            if (k == RM) return e(T(0), -inf<T>, UR);
            if (k == HM) return cinf(-std::sin(y), -std::cos(y));
            if (k == DM) return NN<T>();
            return zz<T>();));
    }
    for (T nu : {T(2), T(3), T(0.3), T(0.5)})
        G(("cyl_bessel_i nu=" + std::to_string((double)nu)).c_str(), quat::cyl_bessel_i(nu, q), TB(isint(nu),
            if (k == FIN) return isint(nu) && y == 0 && x == 0 ? e((std::signbit(x) && int(nu) % 2) ? -T(0) : T(0), y) : skip<T>();
            if (k == RP) return e(inf<T>, y);
            if (k == RM) { const R<T> r = cinf(cospi_(nu), sinpi_(nu)); return sinpi_(nu) == 0 ? e(r.re, y) : r; }
            if (k == HP) return cinf(std::cos(y), std::sin(y));
            if (k == HM) return cinf(cospi_(nu) * std::cos(y) + sinpi_(nu) * std::sin(y), sinpi_(nu) * std::cos(y) - cospi_(nu) * std::sin(y));
            if (k == VV) return zz<T>();
            return NN<T>();));
    for (int kd : {1, 2}) {
        const T nu = T(0.3);
        const auto hk = [&](T x, T y) -> R<T> {
            if (std::isnan(x) || std::isnan(y)) return NN<T>();
            if (fin(x) && fin(y)) return skip<T>();
            if (kd == 1 || !std::isinf(y)) return zz<T>();
            return fin(x) ? cphase(nu / 2, x) : NN<T>();
        };
        G(("cyl_hankel_" + std::to_string(kd)).c_str(), kd == 1 ? quat::cyl_hankel_1(nu, q) : quat::cyl_hankel_2(nu, q), hk(x, y));
        G(("sph_hankel_" + std::to_string(kd)).c_str(), kd == 1 ? quat::sph_hankel_1(2u, q) : quat::sph_hankel_2(2u, q),
          (std::isinf(y) && fin(x) && kd == 2) ? cphase(T(1), x) : hk(x, y));
    }
    for (unsigned n : {0u, 1u, 2u}) {
        G(("sph_bessel n=" + std::to_string(n)).c_str(), quat::sph_bessel(n, q), TB(true,
            if (k == FIN) return y == 0 && x == 0 && n > 0 ? e((std::signbit(x) && n % 2) ? -T(0) : T(0), y) : skip<T>();
            if (k == RP || k == RM) return e(T(0), y, UR);
            if (k == HP || k == HM) return zz<T>();
            if (k == VV) return cphase(T(n) / 2, x);
            return NN<T>();));
        G(("sph_neumann n=" + std::to_string(n)).c_str(), quat::sph_neumann(n, q), TB(true,
            if (k == FIN) return y == 0 && x == 0 ? e(std::signbit(x) && n % 2 == 0 ? inf<T> : -inf<T>, y) : skip<T>();
            if (k == RP || k == RM) return e(T(0), y, UR);
            if (k == HP || k == HM) return zz<T>();
            if (k == VV) return cphase(T(n + 1) / 2, x);
            return NN<T>();));
    }
    G("airy_ai", quat::airy_ai(q), TB(true,
        if (k == FIN) return skip<T>();
        if (k == RP || k == RM) return e(T(0), y, UR);
        if (k == HP || k == DP) return zz<T>();
        return NN<T>();));
    G("airy_bi", quat::airy_bi(q), TB(true,
        if (k == FIN) return skip<T>();
        if (k == RP) return e(inf<T>, y);
        if (k == RM) return e(T(0), y, UR);
        return NN<T>();));
    // ---- elliptic integrals
    const auto zero_t = [&](T x, T y, bool realnan) { return TB(realnan,
        if (k == FIN) return skip<T>();
        return y == 0 ? e(T(0), y, UR | UI) : zz<T>();); };
    G("comp_ellint_1", quat::comp_ellint_1(q), zero_t(x, y, false));
    G("comp_ellint_3", quat::comp_ellint_3(T(0.5), q), zero_t(x, y, false));
    G("ellint_3", quat::ellint_3(T(0.5), q, T(1)), zero_t(x, y, false));
    G("comp_ellint_2", quat::comp_ellint_2(q), TB(false,
        if (k == FIN) return skip<T>();
        return e(std::isinf(y) ? inf<T> : y, -x, y == 0 ? UR : 0);));
    for (T kk : {T(0.5), T(0), T(1), T(-0.5)}) {
        const T K = quat::comp_ellint_1(Q(kk)).w(), Kp = quat::comp_ellint_1(Q(std::sqrt((1 - kk) * (1 + kk)))).w();
        G(("ellint_1 k=" + std::to_string((double)kk)).c_str(), quat::ellint_1(kk, q), TB(true,
            if (k == FIN) return skip<T>();
            if (kk == 0 || k == RP || k == RM) return e(x, y);
            if (k == HP || k == HM) return e(x, nan<T>);
            if (std::isinf(x)) return e(x, Kp);
            const T n = std::round(x / pi<T>);
            return n == 0 ? e(T(0), Kp, UR) : e(2 * K * n, Kp);));
        G(("ellint_2 k=" + std::to_string((double)kk)).c_str(), quat::ellint_2(kk, q), TB(true,
            if (k == FIN) return skip<T>();
            if (kk == 0 || k == RP || k == RM) return e(x, y);
            if (k == HP || k == HM) return e(x, nan<T>);
            if (k == VV) { const T s = std::fmod(std::round(x / pi<T>), T(2)) == 0 ? T(1) : T(-1); return cinf(s * std::sin(x), s * std::cos(x)); }
            return NN<T>();));
    }
    // ---- exponential, sine/cosine and Fresnel integrals
    G("expint", quat::expint(q), TB(true,
        if (k == FIN) return skip<T>();
        if (k == RP) return e(inf<T>, y);
        if (k == RM) return e(T(0), y, UR);
        if (k == HP) return cinf(std::cos(y), std::sin(y));
        if (k == DP) return NN<T>();
        return e(T(0), pi<T>, UR);));
    for (unsigned n : {0u, 1u, 2u, 3u})
        G(("expint n=" + std::to_string(n)).c_str(), quat::expint(n, q), TB(false,
            if (k == FIN) return skip<T>();
            if (k == RP) return e(T(0), y, UR);
            if (k == RM) return n == 0 ? e(-inf<T>, y) : e(-inf<T>, n == 1 ? -pi<T> : -inf<T>);
            if (k == HM) return cinf(-std::cos(y), std::sin(y));
            if (k == DM) return NN<T>();
            return zz<T>();));
    G("sinint", quat::sinint(q), TB(true,
        if (k == FIN) return y == 0 && x == 0 ? e(x, y) : skip<T>();
        if (k == RP) return e(h, y);
        if (k == RM) return e(-h, y);
        if (k == HP) return e(h, T(0), UI);
        if (k == HM) return e(-h, T(0), UI);
        if (k == VV) return cinf(std::sin(x), std::cos(x));
        return NN<T>();));
    G("cosint", quat::cosint(q), TB(false,
        if (k == FIN) return y == 0 && x == 0 ? e(-inf<T>, std::signbit(x) ? pi<T> : y) : skip<T>();
        if (k == RP) return e(T(0), y, UR);
        if (k == HP) return zz<T>();
        if (k == RM || k == HM) return e(T(0), pi<T>, UR);
        if (k == VV) return x == 0 ? e(inf<T>, h) : cinf(std::cos(x), -std::sin(x));
        return NN<T>();));
    for (bool s : {true, false})
        G(s ? "fresnel_s" : "fresnel_c", s ? quat::fresnel_s(q) : quat::fresnel_c(q), TB(true,
            if (k == FIN) return y == 0 && x == 0 ? e(x, y) : skip<T>();
            if (k == RP) return e(T(0.5), y);
            if (k == RM) return e(T(-0.5), y);
            if (k == VV) return x == 0 ? e(x, s ? T(-0.5) : T(0.5)) : NN<T>();
            if (k == DP) return s ? e(-inf<T>, inf<T>) : e(inf<T>, inf<T>);
            if (k == DM) return s ? e(inf<T>, inf<T>) : e(-inf<T>, inf<T>);
            return NN<T>();));
    // ---- polylogarithm, Lambert W
    for (T s : {T(2), T(1), T(0.5), T(0), T(-1), T(-0.5), T(3)})
        G(("polylog s=" + std::to_string((double)s)).c_str(), quat::polylog(s, q), TB(false,
            if (k == FIN) return y == 0 && x == 0 ? e(x, y) : skip<T>();
            const bool real = k == RM;
            if (!(s > 0)) return real ? e(s == 0 ? T(-1) : T(0), y, s == 0 ? 0 : UR) : e(s == 0 ? T(-1) : T(0), T(0), UI | (s == 0 ? 0 : UR));
            if (k == RM) return e(-inf<T>, y);
            if (k == HM) return e(-inf<T>, T(0), UI);
            const T th = (k == RP || k == HP) ? -pi<T> : k == VV ? -h : k == DP ? -3 * pi<T> / 4 : -pi<T> / 4;
            return s > 1 ? e(-inf<T>, inf<T>) : s == 1 ? e(-inf<T>, -th) : e(-inf<T>, T(0), UI);));
    for (int kb : {0, -1, 1})
        G(("lambert_w k=" + std::to_string(kb)).c_str(), quat::lambert_w(kb, q), ([&]() -> R<T> {
            if (std::isnan(x) || std::isnan(y)) return y == 0 ? e(nan<T>, nan<T>) : NN<T>();
            if (fin(x) && fin(y)) return y == 0 && x == 0 && kb == 0 ? e(x, y) : skip<T>();
            T th;
            if (std::isinf(x) && std::isinf(y)) th = std::copysign(x > 0 ? pi<T> / 4 : 3 * pi<T> / 4, y);
            else if (std::isinf(y)) th = std::copysign(h, y);
            else if (x > 0) th = 0;
            else th = std::copysign(pi<T>, y);
            const T im = th + 2 * pi<T> * T(kb);
            return im == 0 ? (y == 0 ? e(inf<T>, y) : e(inf<T>, T(0), UI)) : e(inf<T>, im);
        })());
    // ---- incomplete gamma, beta
    for (T a : {T(2), T(2.3), T(0.5), T(-1), T(0), T(-0.5), T(3)}) {
        const T ga = std::tgamma(a);
        const auto ginc = [&](T x, T y, bool upper, bool reg) -> R<T> {
            const bool realnan = isint(a) && (reg || a > 0);  // the library's real-axis predicates at NaN
            return TB(realnan,
                if (k == FIN) return skip<T>();
                if (reg && !(a > 0)) return NN<T>();
                if (!upper && a <= 0 && isint(a)) return NN<T>();
                const T g = reg ? T(1) : ga, top = upper ? T(0) : g, sg = upper ? T(1) : T(-1);
                if (k == RP) return top == 0 ? e(T(0), y, UR) : e(top, y);
                if (k == HP || k == DP) return top == 0 ? zz<T>() : e(top, T(0), UI);
                if (k == VV) return a < 1 ? (top == 0 ? zz<T>() : e(top, T(0), UI)) : NN<T>();
                if (k == DM) return NN<T>();
                if (k == HM) {
                    const T c = cospi_(a - 1) * std::cos(y) + sinpi_(a - 1) * std::sin(y);
                    const T s = sinpi_(a - 1) * std::cos(y) - cospi_(a - 1) * std::sin(y);
                    return cinf(sg * c, sg * s);
                }
                if (isint(a) && a > 0) return e(ci(sg * cospi_(a - 1)), y);
                if (isint(a)) { const T n = -a, s = std::fmod(n, T(2)) == 0 ? T(-1) : T(1); return e(s * inf<T>, s * pi<T> / std::tgamma(n + 1)); }
                const T c = cospi_(a), s = sinpi_(a);
                return c == 0 ? e(upper ? g : T(0), ci(-sg * s), upper ? 0 : UR) : e(ci(-sg * c), ci(-sg * s)););
        };
        const std::string as = " a=" + std::to_string((double)a);
        G(("tgamma_upper" + as).c_str(), quat::tgamma_upper(a, q), ginc(x, y, true, false));
        G(("tgamma_lower" + as).c_str(), quat::tgamma_lower(a, q), ginc(x, y, false, false));
        G(("gamma_p" + as).c_str(), quat::gamma_p(a, q), ginc(x, y, false, true));
        G(("gamma_q" + as).c_str(), quat::gamma_q(a, q), ginc(x, y, true, true));
    }
    for (T b : {T(2.5), T(2), T(-1)})
        G(("beta b=" + std::to_string((double)b)).c_str(), quat::beta(q, b), TB(true,
            if (k == FIN) return skip<T>();
            if (!(b > 0)) return NN<T>();
            if (k == RM) return isint(b) ? e(T(0), y, UR) : e(nan<T>, y);
            return k == RP ? e(T(0), y, UR) : zz<T>();));
    // ibeta with integer a, b is a polynomial; other parameters are not derived (NaN)
    const auto ibeta_p = [](unsigned, unsigned, CL z) { return 6.0L * z * z - 8.0L * z * z * z + 3.0L * z * z * z * z; };
    G("ibeta 2,3", quat::ibeta(T(2), T(3), q), oracle<T>(x, y, false, ibeta_p, 0, 0));
    G("ibetac 2,3", quat::ibetac(T(2), T(3), q), oracle<T>(x, y, false, [](unsigned, unsigned, CL z) { return 1.0L - (6.0L * z * z - 8.0L * z * z * z + 3.0L * z * z * z * z); }, 0, 0));
    G("ibeta 0.5,2", quat::ibeta(T(0.5), T(2), q), (std::isnan(x) || std::isnan(y)) ? NN<T>() : (fin(x) && fin(y)) ? skip<T>() : NN<T>());
    // ---- polynomials (oracle), odd ones give signed zeros at +-0
    const auto oddz = [&](T x, T y, long double c1) -> R<T> { return e((std::signbit(x) != (c1 < 0)) ? -T(0) : T(0), y); };
    for (unsigned n : {0u, 1u, 2u, 3u, 4u, 5u}) {
        const std::string ns = " n=" + std::to_string(n);
        const long double sg = ((n / 2) % 2) ? -1.0L : 1.0L;  // sign of the linear coefficient of odd n
        const auto poly = [&](T x, T y, CL (*p)(unsigned, unsigned, CL), bool odd) {
            if (fin(x) && fin(y) && y == 0 && x == 0 && odd) return oddz(x, y, sg);
            return oracle<T>(x, y, true, p, n, 0);
        };
        G(("legendre" + ns).c_str(), quat::legendre(n, q), poly(x, y, p_legendre, n % 2));
        G(("hermite" + ns).c_str(), quat::hermite(n, q), poly(x, y, p_hermite, n % 2));
        G(("chebyshev_t" + ns).c_str(), quat::chebyshev_t(n, q), poly(x, y, p_cheb_t, n % 2));
        G(("chebyshev_u" + ns).c_str(), quat::chebyshev_u(n, q), poly(x, y, p_cheb_u, n % 2));
    }
    G("laguerre 3", quat::laguerre(3u, q), oracle<T>(x, y, true, p_laguerre, 3, 0));
    G("laguerre 1", quat::laguerre(1u, q), oracle<T>(x, y, true, p_laguerre, 1, 0));
    G("assoc_laguerre 2,1", quat::assoc_laguerre(2u, 1u, q), oracle<T>(x, y, true, p_laguerre, 2, 1));
    G("assoc_legendre 4,2", quat::assoc_legendre(4u, 2u, q), oracle<T>(x, y, true, p_assoc_legendre, 4, 2));
    G("assoc_legendre 3,2", quat::assoc_legendre(3u, 2u, q),
      (fin(x) && fin(y) && y == 0 && x == 0) ? oddz(x, y, 15.0L) : oracle<T>(x, y, true, p_assoc_legendre, 3, 2));
    for (unsigned l : {1u, 3u, 4u}) {
        // odd m: the zero vector part follows the library's convention s = +i sqrt(x^2 - 1)
        const auto al = [&](T x, T y) -> R<T> {
            if (std::isnan(x) || std::isnan(y)) return NN<T>();
            if (fin(x) && fin(y)) {
                if (y == 0 && x == 0 && (l - 1) % 2 == 1) return oddz(x, y, l == 4 ? -7.5L : 0);
                return skip<T>();
            }
            if (y == 0) {
                R<T> r = oracle<T>(x, T(0), false, p_assoc_legendre, l, 1);
                r.re = 0; r.fl |= UR;
                return r;
            }
            return oracle<T>(x, y, false, p_assoc_legendre, l, 1);
        };
        G(("assoc_legendre l=" + std::to_string(l) + ",1").c_str(), quat::assoc_legendre(l, 1u, q), al(x, y));
    }
    for (unsigned l : {0u, 2u, 3u})
        for (unsigned m : {0u, 1u}) {
            if (m > l) continue;
            const T c00 = T(0.5) / std::sqrt(pi<T>);
            G(("sph_legendre " + std::to_string(l) + "," + std::to_string(m)).c_str(), quat::sph_legendre(l, m, q), TB(true,
                if (k == FIN) return y == 0 && x == 0 && m > 0 ? e(std::signbit(x) ? T(0) : -T(0), y) : skip<T>();
                if (l == 0) return y == 0 ? e(c00, y) : e(c00, T(0), UI);
                if (k == RP || k == RM) return e(nan<T>, y);
                if (k != VV) return NN<T>();
                const T sg = m % 2 ? T(-1) : T(1);
                const T c = cospi_(T(m) / 2) * std::cos(T(l) * x) + sinpi_(T(m) / 2) * std::sin(T(l) * x);
                const T s = sinpi_(T(m) / 2) * std::cos(T(l) * x) - cospi_(T(m) / 2) * std::sin(T(l) * x);
                return cinf(sg * (x == 0 ? cospi_(T(m) / 2) : c), sg * (x == 0 ? sinpi_(T(m) / 2) : s));));
        }
#undef G
}

// Inputs that crashed or hung before (riemann_zeta built a vector of ~9e18 entries).
template <class T> void former_crashes() {
    using Q = quaternion<T>;
    const T I = inf<T>;
    for (const Q& q : {Q(1.5, I, 0, 0), Q(0, 0, -I, 0), Q(I, I, 0, 0), Q(-I, 0, I, I), Q(2, 0, 0, I)}) {
        const Q z = quat::riemann_zeta(q);
        CHECK(std::isnan(z.w()) || z.w() == 1);
        const Q hz = quat::hurwitz_zeta(q, T(0.5));
        CHECK(std::isnan(hz.w()) || std::isinf(hz.w()) || std::isnan(hz.x()) || std::isinf(hz.x()));
    }
}
}  // namespace

TEST(special_function_special_values) { FOR_EACH_TYPE(tables); std::printf("  %d table entries x axes checked\n", entries); }
TEST(former_crashes) { FOR_EACH_TYPE(former_crashes); }
TEST_MAIN
