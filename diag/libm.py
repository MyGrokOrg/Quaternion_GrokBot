import sys, mpmath as mp
mp.mp.prec = 200
F = {"exp": mp.exp, "expm1": mp.expm1, "log": mp.log, "log1p": mp.log1p, "sin": mp.sin, "cos": mp.cos,
     "sinh": mp.sinh, "cosh": mp.cosh, "atan2": mp.atan2, "pow": mp.power, "lgamma": lambda x: mp.log(abs(mp.gamma(x))),
     "tgamma": mp.gamma, "sqrt": mp.sqrt, "cbrt": lambda x: mp.sign(x) * mp.cbrt(abs(x)), "atan": mp.atan}
def h(s): return mp.mpf(float.fromhex(s)) if False else mp.mpf(mp.mpmathify(0)) + parse(s)
def parse(s):
    neg = s.startswith("-"); s = s.lstrip("-")
    m, e = s[2:].split("p"); 
    if "." in m: a, b = m.split(".")
    else: a, b = m, ""
    v = mp.mpf(int(a + b, 16)) * mp.mpf(2) ** (int(e) - 4 * len(b))
    return -v if neg else v
worst = {}
eps = mp.mpf(2) ** -63
for line in open(sys.argv[1]):
    p = line.split()
    if len(p) < 3 or p[0] not in F: continue
    args = [parse(a) for a in p[1:-1]]; got = parse(p[-1]); want = F[p[0]](*args)
    err = abs(got - want) / abs(want) / eps if want != 0 else abs(got) / eps
    if err > worst.get(p[0], (0,))[0]: worst[p[0]] = (err, p[1:-1])
for k, (e, a) in sorted(worst.items(), key=lambda kv: -kv[1][0]):
    print("%-7s worst %8.3g ulp(x87 eps) at %s" % (k, e, " ".join(a)))
