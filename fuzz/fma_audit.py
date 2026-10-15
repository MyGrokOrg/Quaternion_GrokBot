#!/usr/bin/env python3
"""FMA-contraction audit: evaluate every function on the same stress inputs with several
builds of fuzz_eval (typically non-FMA and -march=haswell, gcc and clang), compare each
with mpmath, and report per function the worst and mean error (in eps of double), the
fraction of points whose result bits differ between builds, and a verdict.

  python3 fuzz/fma_audit.py --exe gcc=fe_gcc --exe gcc-hsw=fe_gcc_hsw ... [--n 300]
         [--seed 1] [--only fn1,fn2] [--report audit.md] [--jobs 8]

A function is flagged FMA-SENSITIVE when the worst error of the FMA builds exceeds both 8 eps
and twice the worst of the non-FMA builds, or their mean exceeds the non-FMA mean by 30% (or
the same with the roles exchanged: plain builds losing what the contracted fma keeps).
--sites runs the targeted near-equal-subtraction sites (SITES) instead of whole functions.
Builds whose label ends in "-hsw" are the FMA builds.
"""
import argparse, math, multiprocessing as mpc, os, random, resource, signal, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fuzz as F
from mpmath import mp

mp.prec = 160
EPS = float(F.EPS["double"])


def u(rng, a, b): return rng.uniform(a, b)


def near1(rng):  # 1 - 10^-e (complex offset), either sign
    d = 10 ** -u(rng, 2, 14) * complex(math.cos(t := u(rng, 0, math.pi)), math.sin(t))
    return (1 - d) * rng.choice([1, -1])


# extra functions: (p1, p2, z) samplers and references
def ex_regions(fn, rng):
    r = rng.random()
    if fn in ("comp_ellint_1", "comp_ellint_2"):
        z = near1(rng) if r < 0.5 else complex(u(rng, -3, 3), u(rng, 0, 3))
        return 0, 0, z, "k near +-1" if r < 0.5 else "complex"
    if fn == "comp_ellint_3":
        k = u(rng, -0.99, 0.99)
        z = near1(rng) if r < 0.4 else complex(u(rng, -3, 3), u(rng, 0, 3))
        return k, 0, z, "nu near +-1" if r < 0.4 else "complex"
    if fn in ("ellint_1", "ellint_2"):
        k = rng.choice([u(rng, -0.999, 0.999), 1 - 10 ** -u(rng, 2, 10)])
        z = complex(math.pi / 2 + 10 ** -u(rng, 1, 10) * rng.choice([1, -1]), 10 ** -u(rng, 1, 10)) if r < 0.4 \
            else complex(u(rng, -3, 3), u(rng, 0, 2))
        return k, 0, z, "phi near pi/2" if r < 0.4 else "complex"
    if fn in ("sinint", "cosint"):
        if fn == "cosint" and r < 0.3:  # real zero of Ci
            return 0, 0, complex(0.6165054856207162 + 10 ** -u(rng, 1, 12), 10 ** -u(rng, 1, 12)), "near Ci zero"
        m = 10 ** u(rng, -6, 1.7); t = u(rng, 0, math.pi)
        return 0, 0, complex(m * math.cos(t), m * math.sin(t)), "complex"
    if fn in ("fresnel_c", "fresnel_s"):
        m = 10 ** u(rng, -5, 0.8); t = u(rng, 0, math.pi)
        return 0, 0, complex(m * math.cos(t), m * math.sin(t)), "complex"
    if fn in ("gamma_p", "gamma_q", "tgamma_lower"):
        a = 10 ** u(rng, -1, 1.5)
        z = complex(a * (1 + u(rng, -0.2, 0.2)), u(rng, 0, 0.3 * a)) if r < 0.5 else complex(u(rng, 0, 3 * a), u(rng, 0, 2 * a))
        return a, 0, z, "transition z~a" if r < 0.5 else "complex"
    if fn == "beta":
        a = u(rng, -3.5, 6)
        if abs(a - round(a)) < 1e-3: a += 0.25
        return a, 0, complex(u(rng, -4, 6), u(rng, 0, 4)), "complex"
    if fn in ("polygamma", "trigamma"):
        n = rng.randint(1, 5) if fn == "polygamma" else 1
        z = complex(u(rng, -30, -0.5), u(rng, 0, 2)) if r < 0.4 else complex(u(rng, -5, 30), u(rng, 0, 30))
        return n, 0, z, "reflection" if r < 0.4 else "complex"
    if fn in ("legendre", "assoc_legendre"):
        l = rng.randint(0, 60); m = rng.randint(0, min(l, 8)) if fn == "assoc_legendre" else 0
        z = complex(near1(rng)) if r < 0.3 else complex(u(rng, -1, 1), 10 ** -u(rng, 0, 8)) if r < 0.6 \
            else complex(u(rng, -3, 3), u(rng, 0, 3))
        return l, m, z, "near +-1" if r < 0.3 else "on (-1,1)" if r < 0.6 else "complex"
    if fn in ("laguerre", "assoc_laguerre"):
        n = rng.randint(0, 50); m = rng.randint(0, 6) if fn == "assoc_laguerre" else 0
        return n, m, complex(u(rng, 0, 4 * n + 4), 10 ** -u(rng, 0, 8) if r < 0.5 else u(rng, 0, 5)), "oscillatory"
    if fn == "hermite":
        n = rng.randint(0, 50); s = math.sqrt(2 * n + 1)
        return n, 0, complex(u(rng, -s, s), 10 ** -u(rng, 0, 8) if r < 0.5 else u(rng, 0, 3)), "oscillatory"
    if fn in ("chebyshev_t", "chebyshev_u"):
        n = rng.randint(0, 100)
        return n, 0, complex(u(rng, -1.2, 1.2), 10 ** -u(rng, 0, 10) if r < 0.6 else u(rng, 0, 1)), "near [-1,1]"
    if fn == "ellint_3":
        nu = complex(u(rng, -3, 0.9), u(rng, 0, 2) * (r < 0.5))
        return u(rng, -0.99, 0.99), u(rng, -1.5, 1.5), nu, "complex"
    if fn == "cyl_hankel_2":
        p1, p2, x, y, tag = F.regions("cyl_hankel_1", rng)
        return p1, 0, complex(x, y), tag
    if fn == "sph_hankel_1":
        m = 10 ** u(rng, -2, 2); t = u(rng, 0, math.pi)
        return rng.randint(0, 20), 0, complex(m * math.cos(t), m * math.sin(t)), "complex"
    raise KeyError(fn)


def ex_reference(fn, p1, p2, z):
    if fn == "comp_ellint_1": return mp.ellipk(z * z)
    if fn == "comp_ellint_2": return mp.ellipe(z * z)
    if fn == "comp_ellint_3":  # Carlson form (mp.ellippi can exhaust memory near nu = 1)
        m = mp.mpf(p1) ** 2
        return mp.elliprf(0, 1 - m, 1) + z / 3 * mp.elliprj(0, 1 - m, 1, 1 - z)
    if fn == "ellint_3":  # Carlson form, |phi| <= pi/2 (p2 = phi, z = nu)
        m, ph = mp.mpf(p1) ** 2, mp.mpf(p2)
        sn, cs = mp.sin(ph), mp.cos(ph)
        return sn * mp.elliprf(cs ** 2, 1 - m * sn ** 2, 1) + z / 3 * sn ** 3 * mp.elliprj(cs ** 2, 1 - m * sn ** 2, 1, 1 - z * sn ** 2)
    if fn == "ellint_1": return mp.ellipf(z, mp.mpf(p1) ** 2)
    if fn == "ellint_2": return mp.ellipe(z, mp.mpf(p1) ** 2)
    if fn == "sinint": return mp.si(z)
    if fn == "cosint": return mp.ci(z)
    if fn == "fresnel_c": return mp.fresnelc(z)
    if fn == "fresnel_s": return mp.fresnels(z)
    if fn == "gamma_p": return mp.gammainc(p1, 0, z, regularized=True)
    if fn == "gamma_q": return mp.gammainc(p1, z, mp.inf, regularized=True)
    if fn == "tgamma_lower": return mp.gammainc(p1, 0, z)
    if fn == "beta": return mp.beta(p1, z)
    if fn == "polygamma": return mp.psi(int(p1), z)
    if fn == "trigamma": return mp.psi(1, z)
    if fn == "legendre": return mp.legendre(int(p1), z)
    if fn == "assoc_legendre":
        l, m = int(p1), int(p2)
        c = mp.taylor(lambda t: mp.legendre(l, t), 0, l)  # exact polynomial coefficients
        d = sum(c[k] * mp.ff(k, m) * z ** (k - m) for k in range(m, l + 1))
        return mp.sqrt(1 - z * z) ** m * d
    if fn == "laguerre": return mp.laguerre(int(p1), 0, z)
    if fn == "assoc_laguerre": return mp.laguerre(int(p1), int(p2), z)
    if fn == "hermite": return mp.hermite(int(p1), z)
    if fn == "chebyshev_t": return mp.chebyt(int(p1), z)
    if fn == "chebyshev_u": return mp.chebyu(int(p1), z)
    if fn == "cyl_hankel_2": return mp.hankel2(p1, z)
    if fn == "sph_hankel_1": return mp.sqrt(mp.pi / (2 * z)) * mp.hankel1(int(p1) + mp.mpf(0.5), z)
    raise KeyError(fn)


EXTRA = ("comp_ellint_1", "comp_ellint_2", "comp_ellint_3", "ellint_1", "ellint_2", "ellint_3", "sinint", "cosint",
         "fresnel_c", "fresnel_s", "gamma_p", "gamma_q", "tgamma_lower", "beta", "polygamma", "trigamma",
         "legendre", "assoc_legendre", "laguerre", "assoc_laguerre", "hermite", "chebyshev_t",
         "chebyshev_u", "cyl_hankel_2", "sph_hankel_1")


def lg(rng, a, b):  # log-uniform magnitude 10^[a, b]
    return 10 ** u(rng, a, b)


def pm(rng): return rng.choice([1, -1])


# Targeted sites: (function, sampler[, cond]) with inputs aimed at one near-equal subtraction
# in the implementation.  Each sampler returns (p1, p2, x, y, tag).  cond = "im" also measures
# the (tiny) imaginary part relative to itself; cond = "deriv", "ell3" (d/d phi) or "airy"
# divides the error by the condition number max(1, |z f'(z) / f(z)|) or max(1, |z|^{3/2})
# (oscillatory recurrences near their zeros, Airy phase), so that the comparison between
# builds is not dominated by the conditioning of the function itself.
SITES = {
    # Hermite integrand (a - it)^{-s} - (a + it)^{-s} over e^{2 pi t} - 1 (fixed in 1adef77)
    "hurwitz Hermite integrand": ("hurwitz_zeta", lambda r: (u(r, 0.1, 5), 0, u(r, -9, 12), u(r, 20, 200) * pm(r), "large |Im s|")),
    # reflection zeta(s) = 2^s pi^{s-1} sin(pi s/2) Gamma(1-s) zeta(1-s) and 1 - 2^{1-s} near s = 1
    "zeta near s = 0, 1 (1 - 2^(1-s))": ("riemann_zeta", lambda r: (0, 0, rng_c(r, 1e-14, 1e-3, r.choice([0.0, 1.0]))[0],
                                                                   rng_c(r, 1e-14, 1e-3, 0)[1], "near 0/1")),
    "zeta reflection (Re s < 0)": ("riemann_zeta", lambda r: (0, 0, u(r, -60, -0.5), u(r, 0, 40) * (r.random() < 0.7), "Re s < 0")),
    # Gamma(z) Gamma(1-z) = pi / sin(pi z) next to the poles
    "tgamma reflection near poles": ("tgamma", lambda r: (0, 0, -r.randint(1, 25) + pm(r) * lg(r, -12, -2), lg(r, -14, -2) * (r.random() < 0.6), "z ~ -n")),
    "lgamma reflection near poles": ("lgamma", lambda r: (0, 0, -r.randint(1, 40) + pm(r) * lg(r, -12, -2), lg(r, -14, -1), "z ~ -n")),
    "lgamma near the zeros 1, 2": ("lgamma", lambda r: (0, 0, r.choice([1.0, 2.0]) + pm(r) * lg(r, -14, -2), lg(r, -14, -2) * (r.random() < 0.6), "z ~ 1, 2")),
    "digamma reflection (pi cot pi z)": ("digamma", lambda r: (0, 0, -r.randint(0, 30) + pm(r) * lg(r, -10, -1), lg(r, -12, -1) * (r.random() < 0.6), "z ~ -n")),
    "digamma near its zero 1.4616": ("digamma", lambda r: (0, 0, 1.4616321449683622 + pm(r) * lg(r, -14, -3), lg(r, -14, -3) * (r.random() < 0.6), "z ~ x0")),
    "polygamma reflection": ("polygamma", lambda r: (r.randint(1, 4), 0, -r.randint(1, 20) + pm(r) * lg(r, -6, -1), lg(r, -8, -1) * (r.random() < 0.6), "z ~ -n")),
    # erf series / continued-fraction split and erfc = 1 - erf
    "erf/erfc region split |Re z| in [0.4, 1.3]": ("erfc", lambda r: (0, 0, u(r, 0.4, 1.3) * pm(r), u(r, 0, 4) * (r.random() < 0.8), "split")),
    "erf near 0 (z - z^3/3)": ("erf", lambda r: (0, 0, *rng_c(r, 1e-12, 1e-2, 0), "tiny z")),
    # Ei(x) near its real zero, E1 series (-gamma - log z + ...)
    "Ei near the zero 0.3725": ("expint", lambda r: (0, 0, 0.37250741078136663 + pm(r) * lg(r, -14, -3), lg(r, -14, -3) * (r.random() < 0.6), "x ~ x0")),
    # W near the branch point -1/e: p = sqrt(2 (e z + 1))
    "lambert_w near -1/e": ("lambert_w", lambda r: (r.choice([0, -1]), 0, -0.36787944117144233 + pm(r) * lg(r, -14, -3), lg(r, -14, -3) * (r.random() < 0.7) or 1e-300, "z ~ -1/e")),
    # Li_s near z = 1 (log z series, zeta(s - k) terms) and near integer s
    "polylog near z = 1": ("polylog", lambda r: (r.choice([float(r.randint(-2, 5)), u(r, 0.5, 5), r.randint(1, 5) + pm(r) * lg(r, -10, -2)]), 0,
                                                 *rng_c(r, 1e-10, 0.3, 1.0), "z ~ 1")),
    # Bessel reflection Y = (J cos nu pi - J_-nu)/sin nu pi, K = pi/2 (I_-nu - I_nu)/sin nu pi
    "cyl_neumann near integer order": ("cyl_neumann", lambda r: (r.randint(0, 6) + pm(r) * lg(r, -12, -2), 0, *rng_c(r, 1e-2, 30, 0), "nu ~ n")),
    "cyl_bessel_k near integer order": ("cyl_bessel_k", lambda r: (r.randint(0, 6) + pm(r) * lg(r, -12, -2), 0, *rng_c(r, 1e-2, 30, 0, right=True), "nu ~ n")),
    "cyl_hankel_1 upper half plane": ("cyl_hankel_1", lambda r: (u(r, -3, 6), 0, u(r, -30, 30), u(r, 2, 40), "Im z > 0")),
    "cyl_bessel_j near real zeros": ("cyl_bessel_j", lambda r: (r.choice([0.0, 1.0, 2.5]), 0, u(r, 0.5, 40), lg(r, -14, -4) * (r.random() < 0.6), "on axis")),
    "sph_bessel small z (cancelling closed form)": ("sph_bessel", lambda r: (r.randint(1, 10), 0, *rng_c(r, 1e-4, 1, 0), "|z| < 1")),
    "airy negative axis": ("airy_ai", lambda r: (0, 0, -u(r, 1, 60), lg(r, -12, -1) * (r.random() < 0.6), "x < 0"), "airy"),
    "airy_bi negative axis": ("airy_bi", lambda r: (0, 0, -u(r, 1, 60), lg(r, -12, -1) * (r.random() < 0.6), "x < 0"), "airy"),
    # incomplete gamma: P = 1 - Q / Q = 1 - P near the transition, Gamma(a, z) recurrence for a <= 0
    "tgamma_upper a <= 0 (recurrence)": ("tgamma_upper", lambda r: (-u(r, 0, 8), 0, *rng_c(r, 1e-2, 5, 0, right=True), "a <= 0")),
    "gamma_p/gamma_q complement": ("gamma_q", lambda r: ((a := lg(r, -1, 1.5)), 0, a * (1 + u(r, -0.05, 0.05)), u(r, 0, 0.05 * a), "z ~ a")),
    # Carlson forms with 1 - k^2 and 1 - nu next to 1
    "comp_ellint_2 k near 1 (1 - k^2)": ("comp_ellint_2", lambda r: (0, 0, 1 - lg(r, -14, -2), lg(r, -14, -3) * (r.random() < 0.5), "k ~ 1")),
    "comp_ellint_1 k near +-1 (1 - k^2)": ("comp_ellint_1", lambda r: (0, 0, (1 - lg(r, -15, -2)) * pm(r), lg(r, -15, -3) * (r.random() < 0.4), "k ~ +-1")),
    "comp_ellint_3 real k near 1 (1 - k^2)": ("comp_ellint_3", lambda r: ((1 - lg(r, -15, -2)) * pm(r), 0, u(r, -3, 0.9), u(r, 0, 2) * (r.random() < 0.5), "k ~ +-1")),
    "ellint_1 k near 1 and phi near pi/2 (1 - k^2 sin^2)": ("ellint_1", lambda r: (1 - lg(r, -15, -2), 0, 1.5707963267948966 + pm(r) * lg(r, -15, -2), lg(r, -15, -3) * (r.random() < 0.4), "k, phi ~ 1, pi/2")),
    "ellint_3 k near 1 and nu sin^2 phi near 1": ("ellint_3", lambda r: ell3_pole(r), "ell3"),
    "ellint_2 phi near pi/2": ("ellint_2", lambda r: (u(r, 0.5, 0.999), 0, 1.5707963267948966 + pm(r) * lg(r, -12, -2), lg(r, -12, -2) * (r.random() < 0.5), "phi ~ pi/2")),
    # Ci(z) = gamma + log z + (cos z - 1)/z... series
    "cosint/sinint small z": ("sinint", lambda r: (0, 0, *rng_c(r, 1e-10, 1e-1, 0), "tiny z")),
    "fresnel_s small z": ("fresnel_s", lambda r: (0, 0, *rng_c(r, 1e-8, 1e-1, 0), "tiny z")),
    # three-term recurrences inside the oscillatory region
    "legendre recurrence near +-1": ("legendre", lambda r: (r.randint(2, 80), 0, (1 - lg(r, -12, -2)) * pm(r), lg(r, -14, -4) * (r.random() < 0.5), "x ~ +-1"), "deriv"),
    "chebyshev_t recurrence near +-1": ("chebyshev_t", lambda r: (r.randint(2, 150), 0, (1 - lg(r, -12, -2)) * pm(r), lg(r, -14, -4) * (r.random() < 0.5), "x ~ +-1"), "deriv"),
    "laguerre recurrence (oscillatory)": ("laguerre", lambda r: ((n := r.randint(5, 60)), 0, u(r, 0.1, 3.5 * n), lg(r, -14, -4) * (r.random() < 0.5), "osc"), "deriv"),
    "hermite recurrence (oscillatory)": ("hermite", lambda r: ((n := r.randint(5, 60)), 0, u(r, -1, 1) * math.sqrt(2 * n + 1), lg(r, -14, -4) * (r.random() < 0.5), "osc"), "deriv"),
    # axis reduction: Im f(a + i|v|) / |v| * v with tiny |v| (elementary lift)
    "exp/log/sqrt axis reduction, tiny |v|": ("log", lambda r: (0, 0, u(r, 0.01, 30) * pm(r), lg(r, -300, -8), "tiny |v|")),
    "asin/acos sqrt(A^2 - 1) near 1": ("acos", lambda r: (0, 0, (1 + pm(r) * lg(r, -14, -3)) * pm(r), lg(r, -14, -3) * (r.random() < 0.7) or lg(r, -300, -20), "z ~ +-1")),
    "acosh near 1": ("acosh", lambda r: (0, 0, 1 + pm(r) * lg(r, -14, -3), lg(r, -14, -3), "z ~ 1")),
    "atanh near +-1 (log1p differences)": ("atanh", lambda r: (0, 0, (1 - lg(r, -14, -3)) * pm(r), lg(r, -14, -3), "z ~ +-1")),
    # special-function lift: Im f(x + i|v|) / |v| for tiny |v| must stay relatively accurate
    # (measured on Im f alone)
    "tgamma tiny |v| (Im f / |v|)": ("tgamma", lambda r: (0, 0, u(r, *(-8.5, 8)), lg(r, -150, -6), "tiny |v|"), "im"),
    "lgamma tiny |v| (Im f / |v|)": ("lgamma", lambda r: (0, 0, u(r, *(0.1, 20)), lg(r, -150, -6), "tiny |v|"), "im"),
    "digamma tiny |v| (Im f / |v|)": ("digamma", lambda r: (0, 0, u(r, *(-6.5, 20)), lg(r, -150, -6), "tiny |v|"), "im"),
    "erf tiny |v| (Im f / |v|)": ("erf", lambda r: (0, 0, u(r, *(-3, 3)), lg(r, -150, -6), "tiny |v|"), "im"),
    "riemann_zeta tiny |v| (Im f / |v|)": ("riemann_zeta", lambda r: (0, 0, u(r, *(-8, 6)), lg(r, -150, -6), "tiny |v|"), "im"),
    "expint tiny |v| (Im f / |v|)": ("expint", lambda r: (0, 0, u(r, *(0.05, 20)), lg(r, -150, -6), "tiny |v|"), "im"),
    "lambert_w tiny |v| (Im f / |v|)": ("lambert_w", lambda r: (0, 0, u(r, *(0.05, 20)), lg(r, -150, -6), "tiny |v|"), "im"),
    "cyl_bessel_j tiny |v| (Im f / |v|)": ("cyl_bessel_j", lambda r: (u(r, 0, 4), 0, u(r, *(0.1, 30)), lg(r, -150, -6), "tiny |v|"), "im"),
    "airy_ai tiny |v| (Im f / |v|)": ("airy_ai", lambda r: (0, 0, u(r, *(-10, 5)), lg(r, -150, -6), "tiny |v|"), "im"),
    "polylog tiny |v| (Im f / |v|)": ("polylog", lambda r: (u(r, 0.5, 4), 0, u(r, *(-5, 0.9)), lg(r, -40, -6), "tiny |v|"), "im"),
    "tan/tanh Kahan formula, small x": ("tanh", lambda r: (0, 0, lg(r, -12, -2) * pm(r), u(r, 0, 3), "small x")),
}


def ell3_pole(r):
    """k next to 1 and 1 - nu sin^2 phi next to 0 (real nu only below the pole: the cut)."""
    k, ph = 1 - lg(r, -15, -2), 1.5707963267948966 - lg(r, -8, 0)
    y = lg(r, -14, -2) * (r.random() < 0.7)
    d = lg(r, -14, -2) * (pm(r) if y else 1)
    return k, ph, (1 - d) / math.sin(ph) ** 2, y, "nu s^2 ~ 1"


def rng_c(r, lo, hi, centre, right=False):
    """centre + a log-uniform complex offset of size [lo, hi] in the upper half plane
    (right half too if right)."""
    m = math.exp(u(r, math.log(lo), math.log(hi)))
    t = u(r, 0, math.pi / 2 if right else math.pi)
    return centre + m * math.cos(t), m * math.sin(t)


def samples(fn, n, seed):
    rng = random.Random("%s/%d" % (fn, seed))
    out = []
    site = SITES.get(fn)
    for _ in range(n):
        if site:
            p1, p2, x, y, tag = site[1](rng)
        elif fn in EXTRA:
            p1, p2, z, tag = ex_regions(fn, rng)
            x, y = z.real, z.imag
        else:
            p1, p2, x, y, tag = F.regions(fn, rng)
        axis = rng.choice([1, 2, 3, -1, -2, -3]) if y != 0 else 1
        if y < 0: y, axis = -y, -axis
        out.append((float(p1), float(p2), float(x), float(y), axis, tag))
    return out


def ref_scale(fn, p1, p2, x, y):
    if fn in EXTRA:
        with mp.workprec(mp.prec + 64):
            r = +ex_reference(fn, p1, p2, mp.mpc(x, y))
        return r, abs(r)
    r = F.reference(fn, p1, x, y)
    return r, F.envelope(fn, p1, x, y, r)


class Timeout(BaseException):  # not swallowed by mpmath's internal "except Exception"
    pass


def _alarm(*_):
    raise Timeout()


def init_worker():
    # A few mpmath references (huge orders or arguments) can take minutes and gigabytes;
    # cap each worker and skip such points instead of stalling the pool.
    resource.setrlimit(resource.RLIMIT_AS, (2 << 30, 2 << 30))
    signal.signal(signal.SIGALRM, _alarm)


def work(args):
    name, n, seed, exes = args
    S = samples(name, n, seed)
    fn = SITES[name][0] if name in SITES else name
    inp = "".join("%s %s %s %s %s %d\n" % (fn, F.hexf(p1), F.hexf(p2), F.hexf(x), F.hexf(y), ax)
                  for p1, p2, x, y, ax, tag in S)
    outs = {lab: subprocess.run([exe, "double"], input=inp, capture_output=True, text=True).stdout.split("\n")
            for lab, exe in exes}
    st = {lab: [0.0, 0.0, 0, None] for lab, _ in exes}  # worst, sum, count, worst point
    ndiff = 0
    spread = [0.0, None]  # largest |FMA - non-FMA| (same compiler) in eps of the scale
    for i, (p1, p2, x, y, ax, tag) in enumerate(S):
        lines = {lab: outs[lab][i] for lab, _ in exes}
        ndiff += len(set(lines.values())) > 1
        try:
            signal.alarm(20)
            if name in SITES and len(SITES[name]) > 2 and SITES[name][2] == "im" and y < 1e-20:
                # f real on the real axis: f(x + iy) = f(x) + i y f'(x) + O(y^2) (mpmath does
                # not resolve such a tiny imaginary part reliably even with extra precision)
                fr = lambda t: mp.re(ref_scale(fn, p1, p2, t, 0)[0])
                ref = mp.mpc(fr(mp.mpf(x)), y * mp.diff(fr, mp.mpf(x)))
                sc = abs(ref)
            else:
                ref, sc = ref_scale(fn, p1, p2, x, y)
        except (Exception, Timeout):
            continue
        finally:
            signal.alarm(0)
        if not mp.isfinite(abs(ref)) or sc == 0 or not (F.RANGE["double"][0] <= abs(ref) <= F.RANGE["double"][1]):
            continue
        ck = SITES[name][2] if name in SITES and len(SITES[name]) > 2 else None
        imref = mp.im(ref) if ck == "im" else None
        if ck == "airy":
            sc *= max(1, abs(mp.mpc(x, y)) ** 1.5)
        elif ck == "deriv":
            zz = mp.mpc(x, y)
            sc *= max(1, abs(zz * mp.diff(lambda w: ex_reference(fn, p1, p2, w), zz) / ref))
        elif ck == "ell3":  # Pi(nu, phi, k): |phi dPi/dphi / Pi|, dPi/dphi = 1/((1 - nu s^2) sqrt(1 - k^2 s^2))
            sn2 = mp.sin(mp.mpf(p2)) ** 2
            sc *= max(1, abs(p2 / ((1 - mp.mpc(x, y) * sn2) * mp.sqrt(1 - mp.mpf(p1) ** 2 * sn2)) / ref))
        vals = {}
        for lab, line in lines.items():
            re, im = (F.H(v) for v in line.split())
            vals[lab] = mp.mpc(re, im)
        for lab in vals:
            base = lab[:-4]
            if lab.endswith("-hsw") and base in vals:
                d = vals[lab] - vals[base]
                dv = float(abs(d) / sc / EPS) if mp.isfinite(abs(d)) else (0.0 if vals[lab] == vals[base] else math.inf)
                if not dv <= spread[0]: spread[:] = [dv, (lab, p1, p2, x, y, tag)]
        for lab, line in lines.items():
            re, im = (F.H(v) for v in line.split())
            e = float(abs(mp.mpc(re, im) - ref) / sc / EPS) if mp.isfinite(re) and mp.isfinite(im) else math.inf
            if imref is not None and e != math.inf and imref != 0:
                e = max(e, float(abs(im - imref) / abs(imref) / EPS))
            s = st[lab]
            s[1] += min(e, 1e6); s[2] += 1
            if not e <= s[0]: s[0], s[3] = e, (p1, p2, x, y, tag)
    return name, len(S), ndiff, st, spread


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", action="append", required=True, help="label=path; label ending in -hsw = FMA build")
    ap.add_argument("--n", type=int, default=300)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--only", default="")
    ap.add_argument("--sites", action="store_true", help="targeted near-equal-subtraction sites (SITES) instead of functions")
    ap.add_argument("--report", default="")
    ap.add_argument("--jobs", type=int, default=os.cpu_count())
    a = ap.parse_args()
    exes = [tuple(e.split("=", 1)) for e in a.exe]
    pool_names = list(SITES) if a.sites else list(F.BOUNDS) + list(EXTRA)
    fns = [f for f in pool_names if not a.only or f in a.only.split(",")]
    res = []
    with mpc.Pool(a.jobs, initializer=init_worker, maxtasksperchild=1) as pool:
        jobs = {fn: pool.apply_async(work, ((fn, a.n, a.seed, exes),)) for fn in fns}
        for fn, j in jobs.items():
            try:
                res.append(j.get(timeout=1800))
                print("done:", fn, file=sys.stderr, flush=True)
            except Exception as e:  # worker crashed or hung (mpmath)
                print("skipped:", fn, repr(e), file=sys.stderr, flush=True)
        pool.terminate()
    order = {f: i for i, f in enumerate(fns)}
    res.sort(key=lambda r: order[r[0]])
    labs = [l for l, _ in exes]
    hdr = ("| %s | bits differ | " % ("site (function)" if a.sites else "function") +
           " | ".join("%s worst / mean" % l for l in labs) + " | max FMA-nonFMA | verdict |")
    lines = ["# FMA audit (double, %d stress points per %s, seed %d)" % (a.n, "site" if a.sites else "function", a.seed), "",
             "Errors in eps of double against mpmath (relative to the fuzz envelope); max FMA-nonFMA is the",
             "largest difference between the -march=haswell and the plain build of one compiler, same scale.", "", hdr,
             "|" + "---|" * (len(labs) + 4)]
    flagged = []
    for fn, n, nd, st, spr in res:
        if not any(st[l][2] for l in labs):
            continue
        cells, nf, fm = [], [], []
        for l in labs:
            w, s, c, _ = st[l]
            mean = s / c if c else float("nan")
            cells.append("%.3g / %.2g" % (w, mean))
            (fm if l.endswith("-hsw") else nf).append((w, mean))
        wn = max(w for w, _ in nf); mn = max(m for _, m in nf)
        wf = max(w for w, _ in fm); mf = max(m for _, m in fm)
        # either direction: an FMA build worse than the plain builds, or better (the plain
        # builds lose what fma(a, b, c) keeps, e.g. 1 - k*k)
        sens = (wf > max(8, 2 * wn) or mf > 1.3 * mn + 0.05 or
                wn > max(8, 2 * wf) or mn > 1.3 * mf + 0.05)
        verdict = "FMA-SENSITIVE" if sens else ("insensitive" if nd else "bit-identical")
        if sens: flagged.append(fn)
        name = "%s (%s)" % (fn, SITES[fn][0]) if fn in SITES else fn
        lines.append("| %s | %d/%d | %s | %.3g | %s |" % (name, nd, n, " | ".join(cells), spr[0], verdict))
    for fn, n, nd, st, spr in res:
        if fn in flagged:
            lines.append("")
            lines.append("* %s worst points: " % fn + "; ".join("%s %.3g at %r" % (l, st[l][0], st[l][3]) for l in labs))
    text = "\n".join(lines) + "\n"
    print(text)
    if a.report: open(a.report, "w").write(text)


if __name__ == "__main__":
    main()
