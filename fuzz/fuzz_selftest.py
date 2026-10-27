#!/usr/bin/env python3
"""Checks of fuzz.py's error measure at points the nightly fuzz job has reported.

Run before the fuzz itself (fuzz.yml); exits non-zero if a check fails.
"""
import os
import sys

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fuzz  # noqa: E402

H = float.fromhex
EPS = fuzz.EPS["double"]


def err(fn, p1, x, y, got):
    ref = fuzz.reference(fn, p1, x, y)
    return float(abs(got - ref) / fuzz.envelope(fn, p1, x, y, ref) / EPS), fuzz.bound(fn, p1, x, y, "complex"), ref


def main():
    mp.mp.dps = 40
    bad = 0
    # hurwitz_zeta next to a zero: MSVC x64's result in fuzz run 39261302232 (seed 1123492188),
    # 156 eps relative to |zeta| = 0.00238 (bound 64); relative to L / 16 (L = 1.22, the larger
    # leading term) it is 4.9 eps.
    s, a = H("0x1.ce096f9559e40p-3"), H("0x1.aa1b429d0a901p-2")
    e, b, _ = err("hurwitz_zeta", a, s, 0.0, mp.mpf("0.0023812422150964462908"))
    print("hurwitz_zeta near a zero: %.3g eps (bound %g)" % (e, b))
    bad += not (e <= b and 3 < e < 6)
    # ... but not loosened elsewhere: the 486-eps result of the old kernel (GCC 14, x86-64;
    # clang-cl ARM64 gave 626 eps in run 39261302232, seed 1866713422) at s = 0.164 + 4.600i,
    # a = 1.0076 (|zeta| = 0.64, L / 16 = 0.062) must still fail.
    s, y, a = H("0x1.4fec61e281340p-3"), H("0x1.266bd97b1e700p+2"), H("0x1.01f448092574ap+0")
    e, b, _ = err("hurwitz_zeta", a, s, y, mp.mpc("0.61447936554232629458", "0.18619769681680237161"))
    print("hurwitz_zeta away from zeros, old kernel: %.4g eps (bound %g)" % (e, b))
    bad += not (485 < e < 487 and e > b)
    print("%d failure(s)" % bad)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
