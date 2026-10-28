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
    # airy_bi next to a complex zero of Bi: Apple Clang arm64's result in fuzz run 39122864880
    # (seed 2852262293, n = 250) at z = 5.4274 + 9.6021i, 0.0209 from the zero 5.41394 + 9.58610i
    # (|Bi| = 0.0304 < |Ai| / 4 = 0.054; condition number 526): 210 eps relative to |Bi| (bound
    # 128), 29.5 eps relative to |Ai|.
    x, y = H("0x1.5b5a8e84bebbcp+2"), H("0x1.3344947ca47eep+3")
    got = mp.mpc("0.02854584518403250476", "0.010452693561545449777")
    e, b, ref = err("airy_bi", 0.0, x, y, got)
    print("airy_bi near a zero of Bi: %.4g eps (bound %g)" % (e, b))
    bad += not (e <= b and 29 < e < 30)
    # ... a result 300 eps (of |Ai|) off at the same point must still fail ...
    ai = abs(mp.airyai(mp.mpc(x, y)))
    e, b, _ = err("airy_bi", 0.0, x, y, ref + 300 * EPS * ai)
    print("airy_bi near a zero, 300 eps of |Ai| off: %.4g eps (bound %g)" % (e, b))
    bad += not (e > b)
    # ... and so must a 1e-13 relative error away from the zeros (z = 3 + 2i: |Bi| > |Ai| / 4).
    x, y = 3.0, 2.0
    ref = fuzz.reference("airy_bi", 0.0, x, y)
    e, b, _ = err("airy_bi", 0.0, x, y, ref * (1 + 1e-13))
    print("airy_bi away from zeros, 1e-13 off: %.4g eps (bound %g)" % (e, b))
    bad += not (e > b and fuzz.envelope("airy_bi", 0.0, x, y, ref) == abs(ref))
    print("%d failure(s)" % bad)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
