#!/usr/bin/env python3
"""Randomized accuracy check of the special functions against mpmath.

    fuzz.py --exe build/quaternion_fuzz_eval [--seed S] [--n N] [--types double,"long double"]
            [--only fn1,fn2] [--report report.md]

Inputs are exact binary floating-point numbers (so the comparison measures the library, not
input rounding).  Each function is sampled over several regions: exactly on the real axis,
small |z|, moderate and large |z|, and, for Bessel functions, small, negative, half-integer
and large orders.  The point a + b e is put on a random axis e = +-i, +-j, +-k, which also
exercises the axis reduction.  The error is |got - ref| / scale in units of eps(T), where
scale is |ref| except where the documented accuracy is absolute (oscillatory functions near
their zeros, see ENVELOPE below).  Any error above the per-function bound in BOUNDS is
reported and makes the script exit with status 1.  The seed is printed so a failing run can
be reproduced exactly.
"""
import argparse, math, os, random, subprocess, sys, time
import mpmath as mp

EPS = {"double": 2.0 ** -52, "long double": 2.0 ** -63}
# Results must lie well inside the normal range of the type to be compared.  The upper
# limit leaves 2^12 of headroom: the Bessel kernels form intermediates (H_{nu+1}, K_{nu+1},
# products with e^{i nu pi/2}) up to ~|2 nu / z| larger than the result, so results within a
# few thousand times the largest finite value may overflow (a documented limit).
RANGE = {"double": (mp.mpf(2) ** (-1022 + 64), mp.mpf(2) ** 1012),
         "long double": (mp.mpf(2) ** (-16382 + 64), mp.mpf(2) ** 16372)}

# Bounds in eps, from README.md (accuracy tables) with headroom; (bound, note).
BOUNDS = {
    # elementary functions (quaternion.hpp); README: a few ulp.  The "abs" variants are
    # measured against max(|f|, envelope) where the documented accuracy is absolute.
    "exp": 16, "log": 16, "log10": 16, "sqrt": 16, "pow": 64, "sin": 16, "cos": 16, "tan": 16,
    "sinh": 16, "cosh": 16, "tanh": 16, "asin": 16, "acos": 16, "atan": 16, "asinh": 16,
    "acosh": 16, "atanh": 16,
    "tgamma": 64, "lgamma": 64, "digamma": 64, "erf": 64, "erfc": 64,
    "riemann_zeta": 64, "hurwitz_zeta": None,  # 64 * max(1, |s|): conditioning ~ |s|
    "polylog": 512, "lambert_w": 16, "expint": 64, "expint_n": 64, "tgamma_upper": 128,
    "cyl_bessel_j": 64, "cyl_bessel_j_int": 64, "cyl_neumann": 64, "cyl_bessel_i": 64,
    "cyl_bessel_k": 64, "cyl_hankel_1": 64, "sph_bessel": 64, "sph_neumann": 64,
    "airy_ai": 128, "airy_bi": 128,  # condition ~ |z|^(3/2) up to |z| = 20
}
LARGE_ORDER_BOUND = 2048  # README: up to about 650 ulp for 20 <= nu <= 250


def H(s):
    """Parse a C hex float (any mantissa length) exactly."""
    s = s.strip()
    if s in ("nan", "-nan", "inf", "-inf"):
        return mp.mpf(s.lstrip("-")) * (-1 if s.startswith("-") else 1)
    neg = s.startswith("-")
    s = s.lstrip("+-")
    mant, ex = s[2:].split("p")
    ip, _, fp = mant.partition(".")
    v = mp.mpf(int(ip + fp, 16)) * mp.mpf(2) ** (int(ex) - 4 * len(fp))
    return -v if neg else v


def hexf(x):
    return float(x).hex()


def polar(rng, rmin, rmax, real_axis=False, logr=True, phases=None):
    r = math.exp(rng.uniform(math.log(rmin), math.log(rmax))) if logr else rng.uniform(rmin, rmax)
    if real_axis:
        return (r if rng.random() < 0.5 else -r), 0.0
    t = rng.uniform(-math.pi, math.pi) if phases is None else rng.uniform(*phases)
    return r * math.cos(t), r * math.sin(t)


def box(rng, xr, yr):
    return rng.uniform(*xr), rng.uniform(*yr)


def nu_sample(rng, large=False):
    if large:
        return rng.choice([rng.uniform(20, 100), float(rng.randint(20, 100)), rng.randint(20, 100) + 0.5])
    return rng.choice([rng.uniform(0, 3), rng.uniform(-3, 0), float(rng.randint(0, 12)), rng.randint(0, 6) + 0.5,
                       -(rng.randint(0, 5) + 0.5), rng.uniform(3, 12)])


def regions(fn, rng):
    """Yield (p1, p2, x, y, tag) samples for fn."""
    ra = rng.random() < 0.25  # exactly on the real axis
    tag = "real-axis" if ra else "complex"
    if fn in ELEMENTARY:
        kind = rng.random()
        if kind < 0.2:  # tiny vector part: Im f / |v| in the axis reduction
            x, _ = polar(rng, 1e-3, 30, True)
            y = math.exp(rng.uniform(math.log(1e-300), math.log(1e-6)))
            tag = "tiny |v|"
        elif kind < 0.35:  # next to the branch points / 1 / zeros
            c = rng.choice([1.0, -1.0, 0.0, math.pi, math.pi / 2])
            x, y = polar(rng, 1e-12, 1e-2, ra)
            x += c
            tag += " near %.3g" % c
        else:
            x, y = polar(rng, 1e-4, 700 if fn in ("exp", "sinh", "cosh") else 1e6, ra)
        if fn in ("exp", "sinh", "cosh", "sin", "cos", "tan", "tanh") and abs(x) > 700:
            x = math.copysign(700, x)
        if fn in ("sin", "cos", "tan") and abs(y) > 700: y = math.copysign(700, y)
        if ra:  # stay off the branch cuts on the real axis
            if fn in ("log", "log10", "sqrt", "pow"): x = abs(x)
            if fn in ("asin", "acos", "atanh") and abs(x) > 1: x = math.copysign(rng.uniform(0, 1), x)
            if fn == "acosh" and x < 1: x = 1 + abs(x)
        return (rng.uniform(-4, 4) if fn == "pow" else 0), 0, x, y, tag
    if fn == "tgamma":
        x, y = box(rng, (-30, 30), (0, 0) if ra else (-30, 30)); return 0, 0, x, y, tag
    if fn in ("lgamma", "digamma"):
        x, y = box(rng, (0.01, 60) if ra else (-60, 60), (0, 0) if ra else (-60, 60)); return 0, 0, x, y, tag
    if fn in ("erf", "erfc"):
        x, y = polar(rng, 1e-4, 6, ra); return 0, 0, x, y, tag
    if fn == "riemann_zeta":
        if rng.random() < 0.2:  # next to s = 0 and s = 1
            x, y = polar(rng, 1e-12, 0.3, ra)
            if rng.random() < 0.5: x += 1.0
            return 0, 0, x, y, tag + " near 0/1"
        x, y = box(rng, (-2, 3), (0, 0) if ra else (-300, 300)); return 0, 0, x, y, tag
    if fn == "hurwitz_zeta":
        x, y = box(rng, (-10, 12), (0, 0) if ra else (-60, 60)); return rng.uniform(0.05, 5), 0, x, y, tag
    if fn == "polylog":
        s = rng.choice([float(rng.randint(-3, 6)), rng.uniform(0.5, 6),
                        rng.randint(1, 6) + rng.choice([-1, 1]) * math.exp(rng.uniform(math.log(1e-10), math.log(0.2)))])
        x, y = polar(rng, 1e-3, 6, ra)
        if ra and x > 1: x = -x  # real cut (1, inf)
        return s, 0, x, y, tag
    if fn == "lambert_w":
        k = rng.choice([0, 0, -1, 1])
        x, y = polar(rng, 0.5, 50, ra and k == 0)
        if ra and k == 0 and x < 0: x = -x  # branch cut (-inf, -1/e]
        return k, 0, x, y, tag
    if fn == "expint":
        x, y = polar(rng, 1e-3, 40, ra)
        if ra and x < 0: x = -x  # cut on the negative axis
        return 0, 0, x, y, tag
    if fn == "expint_n":
        x, y = polar(rng, 1e-2, 40, ra, phases=(-1.5, 1.5))
        if ra: x = abs(x)
        return rng.randint(1, 6), 0, x, y, tag
    if fn == "tgamma_upper":
        x, y = polar(rng, 1e-2, 100, ra, phases=(-3.0, 3.0))
        if ra: x = abs(x)
        return rng.uniform(-8, 20), 0, x, y, tag
    if fn in ("cyl_bessel_j", "cyl_neumann", "cyl_bessel_i", "cyl_bessel_k", "cyl_hankel_1"):
        large = rng.random() < 0.15
        nu = nu_sample(rng, large)
        lo, hi = rng.choice([(1e-3, 2), (2, 16), (16, 60), (60, 400)])
        x, y = polar(rng, lo, hi, ra)
        return nu, 0, x, y, tag + (" large-order" if large else "") + " |z|~%g" % hi
    if fn == "cyl_bessel_j_int":
        n = rng.choice([rng.randint(-12, 12), rng.randint(20, 80)])
        lo, hi = rng.choice([(1e-3, 2), (2, 16), (16, 60), (60, 400)])
        x, y = polar(rng, lo, hi, ra)
        return n, 0, x, y, tag + (" large-order" if abs(n) >= 20 else "") + " |z|~%g" % hi
    if fn in ("sph_bessel", "sph_neumann"):
        lo, hi = rng.choice([(1e-3, 2), (2, 16), (16, 60)])
        x, y = polar(rng, lo, hi, ra)
        return rng.randint(0, 10), 0, x, y, tag + " |z|~%g" % hi
    if fn in ("airy_ai", "airy_bi"):
        x, y = polar(rng, 1e-3, 20, ra)
        return 0, 0, x, y, tag
    raise KeyError(fn)


# Known accuracy limits that do not depend on FMA contraction (FMA audit): reported in
# the summary, not failed.  log/log10 near |z| = 1 compute log|z| as log(hypot) and are
# absolutely (not relatively) accurate there (log|z| = log(hypot(x, y)) cancels; see README).
KNOWN = {("log", "complex near 1"), ("log10", "complex near 1")}

ELEMENTARY = ("exp", "log", "log10", "sqrt", "pow", "sin", "cos", "tan", "sinh", "cosh", "tanh",
              "asin", "acos", "atan", "asinh", "acosh", "atanh")


def reference(fn, p1, x, y):
    z = mp.mpc(x, y) if y != 0 else mp.mpf(x)
    if fn in ELEMENTARY:
        zc = mp.mpc(x, y)  # complex functions on the upper side (y >= 0) of every cut
        if fn == "pow": return mp.power(zc, p1)
        f = {"log10": mp.log10}.get(fn) or getattr(mp, fn)
        # mpmath loses a tiny imaginary part to cancellation (e.g. atan's log formula):
        # carry enough extra bits to resolve it.
        extra = 64 + (max(0, -int(mp.floor(mp.log(abs(y), 2)))) if y != 0 else 0)
        with mp.workprec(mp.mp.prec + extra):
            return +f(zc)
    if fn == "tgamma": return mp.gamma(z)
    if fn == "lgamma": return mp.loggamma(z)
    if fn == "digamma": return mp.digamma(z)
    if fn == "erf": return mp.erf(z)
    if fn == "erfc": return mp.erfc(z)
    if fn == "riemann_zeta": return mp.zeta(z)
    if fn == "hurwitz_zeta": return mp.zeta(z, p1)
    if fn == "polylog": return mp.polylog(p1, z)
    if fn == "lambert_w": return mp.lambertw(z, int(p1))
    if fn == "expint": return mp.ei(z)
    if fn == "expint_n": return mp.expint(int(p1), z)
    if fn == "tgamma_upper": return mp.gammainc(p1, z)
    if fn in ("cyl_bessel_j", "cyl_bessel_j_int"): return mp.besselj(p1, z)
    if fn == "cyl_neumann": return mp.bessely(p1, z)
    if fn == "cyl_bessel_i": return mp.besseli(p1, z)
    if fn == "cyl_bessel_k": return mp.besselk(p1, z)
    if fn == "cyl_hankel_1":
        # mpmath forms J + iY, which cancels badly in the upper half plane (|H1| ~ e^-Im z
        # while |J| ~ e^Im z); add the lost digits and recompute.
        h = mp.hankel1(p1, z)
        lost = int(mp.log10(abs(mp.besselj(p1, z)) / abs(h) + 1)) if h != 0 else 0
        if lost > 0:
            with mp.extradps(lost + 10):
                h = mp.hankel1(p1, z)
        return h
    if fn in ("sph_bessel", "sph_neumann") and y == 0 and x < 0:
        # j_n and y_n are single-valued: j_n(-x) = (-1)^n j_n(x), y_n(-x) = (-1)^{n+1} y_n(x)
        sgn = (-1) ** (int(p1) + (fn == "sph_neumann"))
        return sgn * reference(fn, p1, -x, y)
    if fn == "sph_bessel": return mp.sqrt(mp.pi / (2 * z)) * mp.besselj(p1 + 0.5, z) if z != 0 else mp.mpf(p1 == 0)
    if fn == "sph_neumann": return mp.sqrt(mp.pi / (2 * z)) * mp.bessely(p1 + 0.5, z)
    if fn == "airy_ai": return mp.airyai(z)
    if fn == "airy_bi": return mp.airybi(z)
    raise KeyError(fn)


def envelope(fn, p1, x, y, ref):
    """Scale for the error: |ref|, or the documented absolute-accuracy envelope."""
    z = mp.mpc(x, y)
    az = abs(z)
    s = abs(ref)
    if az == 0:
        return s
    # Oscillatory J, Y (and j_n, y_n): near zeros the accuracy is relative to the envelope
    # |H| ~ sqrt(2/(pi|z|)) e^{|Im z|} (J for |z| > 16, Y everywhere; README).
    if fn in ("cyl_bessel_j", "cyl_bessel_j_int", "cyl_neumann", "sph_bessel", "sph_neumann"):
        if az > 16 or fn in ("cyl_neumann", "sph_neumann"):
            env = abs(mp.hankel1(p1 + (0.5 if fn.startswith("sph") else 0), z))
            if fn.startswith("sph"):
                env *= mp.sqrt(mp.pi / (2 * az))
            s = max(s, env)
    # I and K: zeros of I on/near the imaginary axis (I(iy) ~ J(y)), and of K left of it.
    if fn in ("cyl_bessel_i", "cyl_bessel_k"):
        if az > 16 or fn == "cyl_bessel_k" or abs(x) < abs(y):
            s = max(s, abs(mp.besselk(p1, z)) if x >= 0 else abs(mp.besseli(p1, z)))
            s = max(s, abs(mp.hankel1(p1, mp.mpc(-y, x))) * (mp.exp(abs(x)) if x < 0 else 1))
    if fn in ("airy_ai", "airy_bi") and x < 0:
        s = max(s, abs(mp.airyai(z)), abs(mp.airybi(z)))
    if fn in ("lgamma",):  # absolute near the zeros of log Gamma (z = 1, 2)
        s = max(s, mp.mpf(1))
    if fn == "expint":  # absolute next to the real zero of Ei (x = 0.3725...)
        s = max(s, min(mp.mpf(1), abs(mp.exp(z)) / 8))
    if fn == "digamma":
        s = max(s, mp.mpf(1))
    if fn == "hurwitz_zeta" and x < 0:
        # Left of 0 (Re s < 0), zeta(s, a) has real zeros where the growing terms
        # (a + n)^-s cancel; the accuracy there is relative to the largest such term.
        n_max = int(mp.ceil(-x)) + 1
        s = max(s, max(abs(mp.power(p1 + n, -z)) for n in range(n_max + 1)))
    if fn == "erf":  # absolute next to the zeros of erf off the axes
        s = max(s, min(abs(mp.erfc(z)), 1))
    if fn == "erfc":
        s = max(s, min(abs(mp.erf(z)) * abs(mp.erfc(z)), abs(mp.exp(-z * z)) / max(az, 1)))
    return s


def bound(fn, p1, x, y, tag):
    if "large-order" in tag:
        return LARGE_ORDER_BOUND
    if fn == "hurwitz_zeta":
        return 64 * max(1.0, abs(complex(x, y)))
    return BOUNDS[fn]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--n", type=int, default=150, help="samples per function and type")
    ap.add_argument("--types", default="double,long double")
    ap.add_argument("--only", default="")
    ap.add_argument("--report", default="")
    a = ap.parse_args()
    seed = a.seed if a.seed is not None else int.from_bytes(os.urandom(4), "little")
    print("fuzz seed: %d   (reproduce with --seed %d --n %d)" % (seed, seed, a.n), flush=True)
    mp.mp.dps = 40
    fns = [f for f in BOUNDS if not a.only or f in a.only.split(",")]
    rng = random.Random(seed)
    samples = []
    for fn in fns:
        for _ in range(a.n):
            p1, p2, x, y, tag = regions(fn, rng)
            axis = rng.choice([1, 2, 3, -1, -2, -3]) if y != 0 else 1
            if y < 0:  # the lift needs b >= 0 on the chosen axis; reverse the axis instead
                y, axis = -y, -axis
            if y != 0 and rng.random() < 0.25:
                # oblique axis (2, 3, 6)/7: vector part (2c, 3c, 6c) with |v| = 7c exactly
                m, e = math.frexp(y / 7)
                c = math.ldexp(round(m * 2 ** 40) / 2 ** 40, e)
                if c > 0 and 7 * c / 7 == c:
                    y, axis = 7 * c, 4 if axis > 0 else -4
            samples.append((fn, float(p1), float(p2), float(x), float(y), axis, tag))
    rows, failures, known = [], [], []
    t0 = time.time()
    refs = {}
    for t in a.types.split(","):
        inp = "".join("%s %s %s %s %s %d\n" % (fn, hexf(p1), hexf(p2), hexf(x), hexf(y), axis)
                      for fn, p1, p2, x, y, axis, tag in samples)
        out = subprocess.run([a.exe, t], input=inp, capture_output=True, text=True, check=True).stdout.split("\n")
        eps = EPS[t]
        worst = {}
        for (fn, p1, p2, x, y, axis, tag), line in zip(samples, out):
            key = (fn, p1, x, y)
            if key not in refs:
                try:
                    r = reference(fn, p1, x, y)
                    refs[key] = (r, envelope(fn, p1, x, y, r))
                except Exception as e:  # mpmath failure: skip the point
                    refs[key] = (None, None)
            ref, sc = refs[key]
            if ref is None or not mp.isfinite(abs(ref)) or sc == 0:
                continue
            if not (RANGE[t][0] <= abs(ref) <= RANGE[t][1]):
                continue  # under/overflows in T
            re, im = (H(v) for v in line.split())
            got = mp.mpc(re, im)
            err = float(abs(got - ref) / sc / eps) if mp.isfinite(re) and mp.isfinite(im) else float("inf")
            if "tiny |v|" in tag and err != float("inf") and 0 < abs(mp.im(ref)) < abs(ref) * mp.mpf(1e-6):
                # the axis reduction must keep a vector part that is tiny because v is tiny
                # (f(x + iy) ~ f(x) + i f'(x) y) relatively accurate
                # (relative down to the smallest normal of T; subnormals are absolute)
                err = max(err, float(abs(im - mp.im(ref)) / max(abs(mp.im(ref)), RANGE[t][0]) / eps))
            b = bound(fn, p1, x, y, tag)
            w = worst.get(fn)
            if w is None or err > w[0]:
                worst[fn] = (err, p1, x, y, tag)
            if not err <= b and (fn, tag) in KNOWN:
                known.append((t, fn, err, tag))
            elif not err <= b:
                failures.append((t, fn, err, b, p1, x, y, tag, re, im, ref))
        for fn in fns:
            if fn in worst:
                err, p1, x, y, tag = worst[fn]
                rows.append((t, fn, err, bound(fn, p1, x, y, tag), p1, x, y, tag))
    lines = ["# mpmath fuzz report", "", "seed `%d`, %d samples per function and type, %.0f s" % (seed, a.n, time.time() - t0), "",
             "| type | function | worst (eps) | bound | at p, z |", "|---|---|---|---|---|"]
    for t, fn, err, b, p1, x, y, tag in rows:
        lines.append("| %s | %s | %.3g | %g | p=%r, z=%r%s%ri (%s) |" % (t, fn, err, b, p1, x, "+" if float(y) >= 0 else "", y, tag))
    if known:
        lines += ["", "## %d out-of-bound point(s) at known non-FMA limits (not failed)" % len(known), ""]
        kw = {}
        for t, fn, err, tag in known:
            kw[(t, fn, tag)] = max(kw.get((t, fn, tag), 0), err)
        lines += ["* %s %s [%s]: worst %.3g eps" % (t, fn, tag, e) for (t, fn, tag), e in sorted(kw.items())]
    lines += ["", "## %d failure(s)" % len(failures), ""]
    for t, fn, err, b, p1, x, y, tag, re, im, ref in failures[:200]:
        lines.append("* %s %s(p=%s, z=%s%s%si) [%s]: %.3g eps > %g; got %s%s%si, mpmath %s" % (
            t, fn, float(p1).hex(), float(x).hex(), "+" if float(y) >= 0 else "", float(y).hex(), tag, err, b,
            mp.nstr(re, 20), "+" if im >= 0 else "", mp.nstr(im, 20), mp.nstr(ref, 20)))
    text = "\n".join(lines) + "\n"
    print(text)
    if a.report:
        open(a.report, "w").write(text)
    if os.environ.get("GITHUB_STEP_SUMMARY"):
        open(os.environ["GITHUB_STEP_SUMMARY"], "a").write(text)
    print("fuzz seed: %d" % seed)
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
