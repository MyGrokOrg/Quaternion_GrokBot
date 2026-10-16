# Generates the constants of include/quaternion/detail/x87_libm.hpp (mpmath, 400 bits):
#   python3 tests/gen/gen_x87_constants.py > /tmp/x87_constants.txt
# Each long double is printed as an exact hex literal 0x<64-bit significand>p<exp>L.
import mpmath as mp
mp.mp.prec = 400

def rn(v, bits=64):
    """v rounded to nearest with a `bits`-bit significand, as an mpf."""
    if v == 0: return mp.mpf(0)
    with mp.workprec(bits): return +v

def lit(v):
    v = mp.mpf(v)
    if v == 0: return "0.0L"
    s = "-" if v < 0 else ""
    v = abs(v)
    e = int(mp.floor(mp.log(v, 2))) - 63
    m = v / mp.mpf(2) ** e
    while m >= 2**64: e += 1; m /= 2
    while m < 2**63: e -= 1; m *= 2
    assert m == int(m), v
    return "%s0x%016xp%dL" % (s, int(m), e)

def dd(name, v):
    hi = rn(v); lo = rn(v - hi)
    print("inline constexpr dd %s{%s, %s};" % (name, lit(hi), lit(lo)))

pi = mp.pi
dd("k_pio2", pi / 2)
dd("k_pi", pi)
dd("k_logpi", mp.log(pi))
dd("k_hl2pi", mp.log(2 * pi) / 2)
dd("k_egamma", mp.euler)
dd("k_1mgamma", 1 - mp.euler)
dd("k_zeta2h", mp.zeta(2) / 2)            # lgamma(1+z): coefficient of z^2
dd("k_zeta2m1h", (mp.zeta(2) - 1) / 2)    # lgamma(2+z): coefficient of z^2
# ln 2 in three parts: L1 with 48 significant bits (k*L1 exact for |k| < 2^16), L2, L3
ln2 = mp.log(2)
L1 = rn(ln2, 48); L2 = rn(ln2 - L1); L3 = rn(ln2 - L1 - L2)
print("inline constexpr ld k_ln2_1 = %s, k_ln2_2 = %s, k_ln2_3 = %s;" % (lit(L1), lit(L2), lit(L3)))
print("inline constexpr ld k_inv_ln2 = %s;" % lit(rn(1 / ln2)))
print("inline constexpr ld k_pio4 = %s;" % lit(rn(pi / 4)))

def arr(name, vals):
    print("inline constexpr ld %s[] = {" % name)
    for v in vals: print("    %s," % lit(rn(v)))
    print("};")
arr("k_sin", [(-1) ** k / mp.factorial(2 * k + 1) for k in range(1, 11)])     # r^3 .. r^21
arr("k_sin_tail", [(-1) ** k / mp.factorial(2 * k + 1) for k in range(2, 11)])
arr("k_cos", [(-1) ** k / mp.factorial(2 * k) for k in range(2, 12)])         # r^4 .. r^22
arr("k_sinh", [1 / mp.factorial(2 * k + 1) for k in range(1, 11)])            # x^3 .. x^21
arr("k_sinh_tail", [1 / mp.factorial(2 * k + 1) for k in range(2, 11)])
arr("k_exp", [1 / mp.factorial(n) for n in range(2, 19)])                     # r^2 .. r^18
arr("k_exp_tail", [1 / mp.factorial(n) for n in range(3, 19)])
arr("k_log_tail", [mp.mpf(1) / (2 * j + 1) for j in range(4, 18)])            # s^4/9 ..: 1/9 .. 1/35
arr("k_lg1", [(-1) ** k * mp.zeta(k) / k for k in range(3, 34)])              # lgamma(1+z): z^3 .. z^33
arr("k_lg2", [(-1) ** k * (mp.zeta(k) - 1) / k for k in range(3, 26)])        # lgamma(2+z): z^3 .. z^25
arr("k_stirling", [mp.bernoulli(2 * k) / (2 * k * (2 * k - 1)) for k in range(1, 12)])  # 1/x .. 1/x^21
arr("k_stirling_tail", [mp.bernoulli(2 * k) / (2 * k * (2 * k - 1)) for k in range(2, 12)])
# 2/pi: bits 1 .. 64*NW after the binary point, 64 per word, most significant first
NW = 262
with mp.workprec(64 * NW + 128):
    t = int(mp.floor(2 / mp.pi * mp.mpf(2) ** (64 * NW)))
words = [(t >> (64 * (NW - 1 - i))) & (2**64 - 1) for i in range(NW)]
print("inline constexpr std::uint64_t k_two_over_pi[%d] = {" % NW)
for i in range(0, NW, 4):
    print("    " + ", ".join("0x%016xu" % w for w in words[i:i + 4]) + ",")
print("};")
