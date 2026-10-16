# Checks fuzz/libm_check_full.cpp output against mpmath (see that file); prints a Markdown table.
# Error in x87 ulps as in fuzz/libm_check.py (|got - exact| / |exact| / 2^-63), plus the share
# of results that are not the correctly rounded 64-bit-significand value.
import sys, mpmath as mp
mp.mp.prec = 256
F = {"exp": mp.exp, "exp2": lambda x: mp.power(2, x), "expm1": mp.expm1, "log": mp.log, "log1p": mp.log1p,
     "sin": mp.sin, "cos": mp.cos, "tan": mp.tan, "sinh": mp.sinh, "cosh": mp.cosh, "tanh": mp.tanh,
     "asin": mp.asin, "acos": mp.acos, "atan": mp.atan, "atan2": mp.atan2, "pow": mp.power,
     "lgamma": lambda x: mp.log(abs(mp.gamma(x))), "tgamma": mp.gamma, "cbrt": lambda x: mp.sign(x) * mp.cbrt(abs(x)),
     "sqrt": mp.sqrt, "hypot": lambda a, b: mp.sqrt(a * a + b * b), "fma": lambda a, b, c: a * b + c}
def parse(s):
    if s in ("inf", "-inf", "nan", "-nan"): return None
    neg = s.startswith("-"); s = s.lstrip("-")
    m, e = s[2:].split("p")
    a, b = m.split(".") if "." in m else (m, "")
    v = mp.mpf(int(a + b, 16)) * mp.mpf(2) ** (int(e) - 4 * len(b))
    return -v if neg else v
def cr64(v):  # round to 64-bit significand, nearest (x87 extended normal range)
    with mp.workprec(64): return +v
eps = mp.mpf(2) ** -63
st = {}
for line in open(sys.argv[1]):
    p = line.split()
    if not p or p[0].startswith("#"): continue
    key = p[0]; f = key.split("@")[0]
    if f not in F: continue
    args = [parse(a) for a in p[1:-1]]; got = parse(p[-1])
    if got is None or any(a is None for a in args): continue
    want = F[f](*args)
    if want == 0 or not mp.isfinite(want) or abs(want) > mp.mpf(2) ** 16383 or abs(want) < mp.mpf(2) ** -16382: continue
    err = abs(got - want) / abs(want) / eps
    s = st.setdefault(key, {"n": 0, "errs": [], "ncr": 0, "worst": (-1, None)})
    s["n"] += 1; s["errs"].append(float(err)); s["ncr"] += got != cr64(want)
    if err > s["worst"][0]: s["worst"] = (err, p[1:-1])
print("| function@range | n | worst (x87 eps) | median | >1 eps | not CR | worst at |")
print("|---|---|---|---|---|---|---|")
for k in sorted(st, key=lambda k: -st[k]["worst"][0]):
    s = st[k]; e = sorted(s["errs"])
    print("| %s | %d | %.3g | %.3g | %d | %d | %s |" % (k, s["n"], s["worst"][0], e[len(e) // 2],
          sum(x > 1 for x in e), s["ncr"], " ".join(s["worst"][1])))
