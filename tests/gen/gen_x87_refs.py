#!/usr/bin/env python3
"""Generate tests/x87_refs.inc: correctly rounded x87 long double references (mpmath) for the
library's own long double kernels (quaternion/detail/x87_libm.hpp, QUATERNION_X87_LIBM).

Each entry is { function, x, y, f(x[, y]) rounded to nearest-even in the 64-bit significand }.
Inputs are exact long doubles drawn with a fixed seed: random significands at every exponent
range, points next to multiples of pi/2 (the hard cases for argument reduction) and special
regions.  Usage: python3 tests/gen/gen_x87_refs.py > tests/x87_refs.inc, and
python3 tests/gen/gen_x87_refs.py --complex > tests/x87_crefs.inc for the std::complex functions
(entries { function, re z, im z, re f(z), im f(z) }, each part correctly rounded; seed + 1).
"""
import random
import sys

import mpmath as mp

mp.mp.prec = 600
SEED = 20261016
N = 160  # per input family
rng = random.Random(SEED)
LDBL_MAX = mp.ldexp((1 << 64) - 1, 16384 - 64)


def to_ld(v):
    """Round an mpf to the nearest long double (normal range, ties to even)."""
    v = mp.mpf(v)
    if v == 0:
        return v
    m, e = mp.frexp(v)  # v = m 2^e, 0.5 <= |m| < 1
    e = max(e, -16381)   # subnormals: fixed exponent
    q = mp.ldexp(v, 64 - e)
    n = int(mp.floor(q))
    f = q - n
    if f > 0.5 or (f == 0.5 and n % 2):
        n += 1
    return mp.ldexp(n, e - 64)


def lit(v):
    v = mp.mpf(v)
    if mp.isinf(v):
        return "inf" if v > 0 else "-inf"
    if v == 0:
        return "0.0L"
    m, e = mp.frexp(v)
    mi = int(mp.ldexp(m, 64))
    s = "-" if mi < 0 else ""
    return "%s0x%016xp%dL" % (s, abs(mi), e - 64)


def rnd_mag(lo, hi, sign=True):
    e = rng.randint(lo, hi - 1)
    m = rng.getrandbits(63) | (1 << 63)
    v = mp.ldexp(m, e - 63)
    return -v if sign and rng.random() < 0.5 else v


def uni(a, b):
    return to_ld(mp.mpf(a) + (mp.mpf(b) - a) * mp.mpf(rng.random()))


def near_pio2(kmax_bits):
    k = rng.getrandbits(kmax_bits) | 1
    return to_ld(k * mp.pi / 2)


FAMILIES = {
    "sin": lambda: [rnd_mag(-70, 2) for _ in range(N)] + [rnd_mag(2, 64) for _ in range(N)]
    + [rnd_mag(64, 16384) for _ in range(N)] + [near_pio2(b) for b in range(1, 121) for _ in range(1)]
    + [near_pio2(rng.randint(1, 30)) * rng.choice([1, -1]) for _ in range(N // 2)],
}
FAMILIES["cos"] = FAMILIES["sin"]
FAMILIES["tan"] = FAMILIES["sin"]


def pow_family():
    out = [(uni(0.01, 100), uni(-50, 50)) for _ in range(N)]
    out += [(mp.mpf(rng.randint(2, 1000)), uni(-8, 8)) for _ in range(N)]
    while len(out) < 3 * N:
        x, y = abs(rnd_mag(-100, 100)), uni(-100, 100)
        if abs(y * mp.log(x)) < 11000:
            out.append((x, y))
    out += [(uni(-10, -0.1), mp.mpf(rng.randint(-30, 30))) for _ in range(N // 2)]
    out += [(uni(0.5, 2), uni(-16000, 16000)) for _ in range(N // 2)]  # results across the range
    return out


def no_poles(xs):
    return [x for x in xs if not (x <= 0 and x == mp.floor(x))]


FAMILIES.update({
    "sinh": lambda: [rnd_mag(-70, 0) for _ in range(N)] + [uni(-50, 50) for _ in range(N)]
    + [uni(-11350, 11350) for _ in range(N)],
    "expm1": lambda: [rnd_mag(-70, 0) for _ in range(N)] + [uni(-50, 50) for _ in range(N)]
    + [uni(-50, 11350) for _ in range(N)],
    "pow": pow_family,
    "tgamma": lambda: no_poles([uni(0.001, 1750) for _ in range(N)] + [uni(-170, 0) for _ in range(N)]
                               + [uni(0.5, 30) for _ in range(N)] + [mp.mpf(n) for n in range(1, 41)]
                               + [abs(rnd_mag(-70, -1)) for _ in range(N // 4)]),
    "lgamma": lambda: no_poles([uni(0.001, 50) for _ in range(N)] + [uni(-50, 0) for _ in range(N)]
                               + [uni(0.7, 2.3) for _ in range(N)] + [abs(rnd_mag(-70, 200)) for _ in range(N)]),
})


def lgamma_ref(x):
    return mp.log(abs(mp.gamma(x))) if x < 0 else mp.loggamma(x)


REF = {"sin": mp.sin, "cos": mp.cos, "tan": mp.tan, "sinh": mp.sinh, "expm1": mp.expm1,
       "tgamma": mp.gamma, "lgamma": lgamma_ref}


CREF = {"exp": mp.exp, "sin": mp.sin, "cos": mp.cos, "sinh": mp.sinh, "cosh": mp.cosh, "sqrt": mp.sqrt}


def complex_family(f):
    out = []
    for _ in range(N):  # moderate: |re| <= 40 for the exponentials, |im| up to 2^20
        a = uni(-40, 40) if f != "sqrt" else rnd_mag(-60, 60)
        b = rnd_mag(-40, 20) if f != "sqrt" else rnd_mag(-60, 60)
        out.append((a, b) if f not in ("sin", "cos") else (b, a))
    for _ in range(N // 4):  # next to multiples of pi/2 in the oscillating part
        a, b = uni(-5, 5), near_pio2(rng.randint(1, 40))
        out.append((a, b) if f not in ("sin", "cos") else (b, a))
    if f == "sqrt":  # next to the negative real axis
        out += [(-abs(rnd_mag(-20, 20)), rnd_mag(-80, -30)) for _ in range(N // 4)]
    return out


def main_complex():
    out = []
    counts = {}
    for f, fn in CREF.items():
        zs = complex_family(f)
        for a, b in zs:
            w = fn(mp.mpc(a, b))
            out.append('    {"%s", %s, %s, %s, %s},' % (f, lit(a), lit(b), lit(to_ld(w.real)), lit(to_ld(w.imag))))
        counts[f] = len(zs)
    sys.stdout.write("// Generated by tests/gen/gen_x87_refs.py --complex (mpmath %s, seed %d, %s); do not edit.\n"
                     % (mp.__version__, SEED + 1, ", ".join("%s %d" % kv for kv in counts.items())))
    sys.stdout.write("\n".join(out) + "\n")


def main():
    out = []
    counts = {}
    for f, fam in FAMILIES.items():
        xs = fam()
        for x in xs:
            x, y = x if isinstance(x, tuple) else (x, mp.mpf(0))
            r = to_ld(mp.power(x, y) if f == "pow" else REF[f](x))
            out.append('    {"%s", %s, %s, %s},' % (f, lit(x), lit(y), lit(r)))
        counts[f] = len(xs)
    sys.stdout.write("// Generated by tests/gen/gen_x87_refs.py (mpmath %s, seed %d, %s); do not edit.\n"
                     % (mp.__version__, SEED, ", ".join("%s %d" % kv for kv in counts.items())))
    sys.stdout.write("\n".join(out) + "\n")


if __name__ == "__main__":
    if sys.argv[1:] == ["--complex"]:
        rng.seed(SEED + 1)
        main_complex()
    else:
        main()
