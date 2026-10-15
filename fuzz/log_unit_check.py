#!/usr/bin/env python3
"""log, log10 and pow next to |q| = 1 against mpmath (README, "log near |z| = 1").

Quaternions q = r u with u a random direction in R^4 (a quarter complex, a quarter with a
dominant real part), |q| = r = 1 +- 10^U(-17, -1); a further quarter has |q| = 10^U(-0.15, 0.15)
and half of those are complex points whose x^2 + y^2 - 1 cancels beyond double length.

  log mode: worst relative error in eps of Re log q = log|q| and Re log10 q (double, long double,
            and float for the inputs rounded to float).
  pow mode: q^p with p uniform in [-4, 4] ("small") or +-10^U(1, 3) ("large"): worst normwise
            relative error, worst relative error of |q^p|, and (with two executables) how many
            results differ in their bits.

Build the evaluator with g++ -std=c++23 -O2 -Iinclude fuzz/log_unit_eval.cpp -o log_unit_eval.
  python3 fuzz/log_unit_check.py log --exe NEW [--exe OLD] --n 20000 --seed 1
  python3 fuzz/log_unit_check.py pow --exe NEW [--exe OLD] --n 20000 --seed 5 --p small
"""
import argparse, math, random, struct, subprocess
from mpmath import mp, mpf, sqrt, log, exp

mp.prec = 300


def H(s):
    if "0x" not in s:
        return mpf(float(s))
    m, e = s.split("p")
    sg = -1 if m.startswith("-") else 1
    ip, _, fp = m.lstrip("-")[2:].partition(".")
    return sg * (mpf(int(ip, 16)) + (mpf(int(fp, 16)) / mpf(16) ** len(fp) if fp else 0)) * mpf(2) ** int(e)


def points(rng, n, mode, pkind):
    pts = []
    for i in range(n):
        kind = i % 4
        v = [rng.gauss(0, 1) for _ in range(4)]
        if kind == 1: v[2] = v[3] = 0.0
        if kind == 2: v[1:] = [x * 10 ** rng.uniform(-9, -1) for x in v[1:]]
        nv = math.sqrt(sum(x * x for x in v))
        r = 1 + rng.choice([-1, 1]) * 10 ** rng.uniform(-17, -1) if kind != 3 else 10 ** rng.uniform(-0.15, 0.15)
        q = [x / nv * r for x in v]
        if kind == 3 and rng.random() < 0.5:
            a = 1 - rng.choice([1, 2, 3, 17]) * 2.0 ** -rng.randint(30, 53)
            q = [a, math.sqrt(1 - a * a) * (1 + rng.uniform(-1e-16, 1e-16)), 0.0, 0.0]
        if mode == "pow":
            q.append(rng.uniform(-4, 4) if pkind == "small" else rng.choice([-1, 1]) * 10 ** rng.uniform(1, 3))
        pts.append(q)
    return pts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["log", "pow"])
    ap.add_argument("--exe", action="append", required=True)
    ap.add_argument("--n", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--p", choices=["small", "large"], default="small")
    a = ap.parse_args()
    pts = points(random.Random(a.seed), a.n, a.mode, a.p)
    inp = "".join(" ".join(float.hex(x) for x in q) + "\n" for q in pts)
    outs = [subprocess.run([e, a.mode], input=inp, capture_output=True, text=True, check=True).stdout.split("\n")[:len(pts)] for e in a.exe]
    eps = [mpf(2) ** -52, mpf(2) ** -63, mpf(2) ** -23]
    print("%s, seed %d, %d points%s" % (a.mode, a.seed, a.n, ", p " + a.p if a.mode == "pow" else ""))
    for e, out in zip(a.exe, outs):
        worst = {}
        for q, line in zip(pts, out):
            f = line.split()
            s = sum(mpf(x) ** 2 for x in q[:4])
            L = log(s) / 2
            if a.mode == "log":
                Lf = log(sum(mpf(struct.unpack("f", struct.pack("f", x))[0]) ** 2 for x in q[:4])) / 2
                for k, (name, ref, t) in enumerate([("log double", L, 0), ("log10 double", L / log(10), 0),
                                                    ("log long double", L, 1), ("log10 long double", L / log(10), 1),
                                                    ("log float", Lf, 2), ("log10 float", Lf / log(10), 2)]):
                    if ref == 0: continue
                    err = float(abs(H(f[k]) - ref) / abs(ref) / eps[t])
                    worst[name] = max(worst.get(name, 0), err)
            else:
                b, c, d = (mpf(x) for x in q[1:4])
                p = mpf(q[4]); vn = sqrt(b * b + c * c + d * d); th = mp.atan2(vn, mpf(q[0]))
                w = exp(p * L)
                ref = [w * mp.cos(p * th)] + [w * mp.sin(p * th) * x / vn for x in (b, c, d)]
                rn = sqrt(sum(x * x for x in ref))
                for t, tn in enumerate(["double", "long double"]):
                    got = [H(x) for x in f[4 * t:4 * t + 4]]
                    err = float(sqrt(sum((g - r) ** 2 for g, r in zip(got, ref))) / rn / eps[t])
                    merr = float(abs(sqrt(sum(g * g for g in got)) - w) / w / eps[t])
                    worst["pow " + tn] = max(worst.get("pow " + tn, 0), err)
                    worst["|pow| " + tn] = max(worst.get("|pow| " + tn, 0), merr)
        print("  %s: " % e + ", ".join("%s %.3g" % kv for kv in worst.items()))
    if len(outs) == 2:
        print("  results that differ in their bits: %d of %d" % (sum(x != y for x, y in zip(*outs)), len(pts)))


if __name__ == "__main__":
    main()
