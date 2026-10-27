// Own long double kernels for the C library functions that Apple's x87 libm (macOS x86_64)
// and mingw-w64's (MinGW-w64 x86_64) get wrong: sinl, cosl, tanl, sinhl, expm1l, powl, tgammal and lgammal (and asinl and acosl
// next to +-1), and the
// std::complex<long double> functions the headers call that libc++ builds on them (exp, sin,
// cos, sinh, cosh, sqrt and polar).  See the README, "Supported toolchains".
//
// Included by quaternion.hpp only when QUATERNION_X87_LIBM is 1: the default on macOS x86_64
// (Apple Clang or LLVM clang) and MinGW-w64 x86_64; -DQUATERNION_X87_LIBM=1 forces it on elsewhere (CI does so on
// Linux x86_64 to test it), -DQUATERNION_X87_LIBM=0 turns it off.  It adds non-template
// overloads for long double to quat::detail::lm, which overload resolution prefers to the
// forwarders there; nothing else changes, and other platforms never see this file.
//
// Method: x87 arithmetic (64-bit significand) and double-length ("dd") values, hi + lo.
//   * sin, cos, tan: Cody-Waite reduction with pi/2 in four parts (232 bits) for
//     pi/4 < |x| < 2^23 pi/2; above, exact Payne-Hanek reduction with 256 bits of 2/pi taken from
//     a 16,768-bit table at the exponent of x (integer arithmetic, so every argument up to
//     LDBL_MAX is reduced to within 2^-190 of a quadrant); then Taylor kernels on |r| <= pi/4.
//   * exp, expm1, sinh: x = k ln2 + r with ln2 in three parts; Taylor kernel for e^r - 1.
//   * pow(x, y) = exp(y log x) with log x double-length (atanh series).
//   * asin, acos for |x| > 1/2: atan2 of a double-length sqrt((1 - |x|)(1 + |x|)).
//   * lgamma: Stirling for x >= 20, upward shift below, Taylor series next to the zeros at 1
//     and 2, reflection for x < 0; tgamma = exp(lgamma) with the sign, exact small factorials.
// Constants: tests/gen/gen_x87_constants.py (mpmath).  Tests: tests/test_x87_libm.cpp.
#ifndef QUATERNION_DETAIL_X87_LIBM_HPP
#define QUATERNION_DETAIL_X87_LIBM_HPP

#if !defined(__x86_64__) || !defined(__GNUC__) || __LDBL_MANT_DIG__ != 64
#error "QUATERNION_X87_LIBM needs the x87 80-bit long double (x86_64, GCC or Clang)"
#endif

#include <cmath>
#include <complex>
#include <cstdint>
#include <bit>
#include <limits>
#include <type_traits>

static_assert(std::numeric_limits<long double>::digits == 64 && std::numeric_limits<long double>::radix == 2 &&
                  std::numeric_limits<long double>::max_exponent == 16384 &&
                  std::numeric_limits<long double>::min_exponent == -16381,
              "QUATERNION_X87_LIBM: long double must be the x87 80-bit format");

namespace quat::detail::lm {
namespace x87 {
using ld = long double;
__extension__ typedef unsigned __int128 u128;

struct dd {
    ld hi, lo;
};
[[nodiscard]] inline dd two_sum(ld a, ld b) noexcept {
    const ld s = a + b, bb = s - a;
    return {s, (a - (s - bb)) + (b - bb)};
}
[[nodiscard]] inline dd fast_two_sum(ld a, ld b) noexcept {  // |a| >= |b| or a == 0
    const ld s = a + b;
    return {s, b - (s - a)};
}
[[nodiscard]] inline dd two_prod(ld a, ld b) noexcept {  // Dekker: x87 has no fma
    constexpr ld c = 4294967297.0L;                       // 2^32 + 1
    const ld p = a * b, ta = c * a, tb = c * b;
    const ld ah = ta - (ta - a), al = a - ah, bh = tb - (tb - b), bl = b - bh;
    return {p, ((ah * bh - p) + ah * bl + al * bh) + al * bl};
}
[[nodiscard]] inline dd neg(dd a) noexcept { return {-a.hi, -a.lo}; }
[[nodiscard]] inline dd add(dd a, dd b) noexcept {
    dd s = two_sum(a.hi, b.hi);
    const dd t = two_sum(a.lo, b.lo);
    s = fast_two_sum(s.hi, s.lo + t.hi);
    return fast_two_sum(s.hi, s.lo + t.lo);
}
[[nodiscard]] inline dd add(dd a, ld b) noexcept {
    const dd s = two_sum(a.hi, b);
    return fast_two_sum(s.hi, s.lo + a.lo);
}
[[nodiscard]] inline dd mul(dd a, ld b) noexcept {
    const dd p = two_prod(a.hi, b);
    return fast_two_sum(p.hi, p.lo + a.lo * b);
}
[[nodiscard]] inline dd mul(dd a, dd b) noexcept {
    const dd p = two_prod(a.hi, b.hi);
    return fast_two_sum(p.hi, p.lo + (a.hi * b.lo + a.lo * b.hi));
}
[[nodiscard]] inline dd div(dd a, dd b) noexcept {
    const ld q1 = a.hi / b.hi;
    dd r = add(a, neg(mul(b, q1)));
    const ld q2 = r.hi / b.hi;
    r = add(r, neg(mul(b, q2)));
    const ld q3 = r.hi / b.hi;
    return add(fast_two_sum(q1, q2), q3);
}
[[nodiscard]] inline ld round(dd a) noexcept { return a.hi + a.lo; }
// x * 2^k for a normal or subnormal result (one rounding: the hardware multiply).
// 2^k for -16382 <= k <= 16383, from the bits (the x87 format: explicit leading significand
// bit, 15-bit biased exponent); cheaper than ldexp, a library call.
struct x87_bits {
    std::uint64_t m;
    std::uint16_t se, pad[3];
};
[[nodiscard]] inline ld pow2(int k) noexcept {
    return std::bit_cast<ld>(x87_bits{std::uint64_t(1) << 63, static_cast<std::uint16_t>(k + 16383), {0, 0, 0}});
}
// biased exponent of a positive normal x (0 for subnormals and zero)
[[nodiscard]] inline int biased_exponent(ld x) noexcept { return std::bit_cast<x87_bits>(x).se & 0x7fff; }
// v rounded to an integer (ties to even) for |v| < 2^62: the x87 significand has 64 bits
[[nodiscard]] inline ld round_int(ld v) noexcept {
    constexpr ld magic = 0x1.8p63L;
    volatile ld t = v + magic;  // volatile: no reassociation even with -fassociative-math
    return t - magic;
}
// x 2^k, rounded once (when the result is subnormal, by the last multiplication)
[[nodiscard]] inline ld scale(ld x, int k) noexcept {
    if (k > 16383) return x * pow2(16383) * pow2(k - 16383);
    if (k >= -16382) return x * pow2(k);
    return x * pow2(k + 256) * pow2(-256);
}
// x 2^k for a double-length x whose parts stay normal (the callers' ranges)
[[nodiscard]] inline dd scale(dd x, int k) noexcept {
    if (k < -16382 || k > 16383) return {std::ldexp(x.hi, k), std::ldexp(x.lo, k)};
    const ld f = pow2(k);
    return {x.hi * f, x.lo * f};
}
template <std::size_t N>
[[nodiscard]] inline ld horner(const ld (&c)[N], ld z) noexcept {
    ld p = c[N - 1];
    for (std::size_t i = N - 1; i-- > 0;) p = p * z + c[i];
    return p;
}

// ---- constants: the output of tests/gen/gen_x87_constants.py, verbatim -------------
inline constexpr dd k_pio2{0xc90fdaa22168c235p-63L, -0xece675d1fc8f8cbbp-129L};
inline constexpr dd k_pi{0xc90fdaa22168c235p-62L, -0xece675d1fc8f8cbbp-128L};
inline constexpr dd k_inv3{0xaaaaaaaaaaaaaaabp-65L, -0xaaaaaaaaaaaaaaabp-130L};
inline constexpr dd k_inv5{0xcccccccccccccccdp-66L, -0xcccccccccccccccdp-132L};
inline constexpr dd k_inv7{0x9249249249249249p-66L, 0x9249249249249249p-132L};
inline constexpr dd k_logpi{0x928682473d0de85fp-63L, -0xa06a93957bc0b668p-128L};
inline constexpr dd k_hl2pi{0xeb3f8e4325f5a535p-64L, -0xd686dffd77cdbfb8p-129L};
inline constexpr dd k_egamma{0x93c467e37db0c7a5p-64L, -0xb90701fbfab4d2a5p-130L};
inline constexpr dd k_1mgamma{0xd8773039049e70b6p-65L, 0xb90701fbfab4d2a5p-130L};
inline constexpr dd k_zeta2h{0xd28d3312983e9918p-64L, 0xe7b1224401759cbdp-129L};
inline constexpr dd k_zeta2m1h{0xa51a6625307d3231p-65L, -0xc276eddff4531a15p-132L};
inline constexpr ld k_ln2_1 = 0xb17217f7d1cf0000p-64L, k_ln2_2 = 0xf35793c7673007e6p-113L, k_ln2_3 = -0x950bf0cbcd98d675p-180L;
inline constexpr ld k_inv_ln2 = 0xb8aa3b295c17f0bcp-63L;
inline constexpr ld k_pio4 = 0xc90fdaa22168c235p-64L;
inline constexpr ld k_2_over_pi = 0xa2f9836e4e44152ap-64L;
inline constexpr ld k_pio2_cw[4] = {0xc90fdaa221000000p-63L, 0xd18469898cc51702p-104L, -0x8f8cbb5bf6c7ddd6p-169L, -0xc19c5fefa20ace27p-234L};
inline constexpr ld k_sin[] = {
    -0xaaaaaaaaaaaaaaabp-66L,
    0x8888888888888889p-70L,
    -0xd00d00d00d00d00dp-76L,
    0xb8ef1d2ab6399c7dp-82L,
    -0xd7322b3faa271c7fp-89L,
    0xb092309d43684be5p-96L,
    -0xd73f9f399dc0f88fp-104L,
    0xca963b81856a5359p-112L,
    -0x97a4da340a0ab926p-120L,
    0xb8dc77b6e7ab8c5fp-129L,
};
inline constexpr ld k_sin_tail[] = {
    0x8888888888888889p-70L,
    -0xd00d00d00d00d00dp-76L,
    0xb8ef1d2ab6399c7dp-82L,
    -0xd7322b3faa271c7fp-89L,
    0xb092309d43684be5p-96L,
    -0xd73f9f399dc0f88fp-104L,
    0xca963b81856a5359p-112L,
    -0x97a4da340a0ab926p-120L,
    0xb8dc77b6e7ab8c5fp-129L,
};
inline constexpr ld k_cos[] = {
    0xaaaaaaaaaaaaaaabp-68L,
    -0xb60b60b60b60b60bp-73L,
    0xd00d00d00d00d00dp-79L,
    -0x93f27dbbc4fae397p-85L,
    0x8f76c77fc6c4bdaap-92L,
    -0xc9cba54603e4e906p-100L,
    0xd73f9f399dc0f88fp-108L,
    -0xb413c31dcbecbbdep-116L,
    0xf2a15d201011283dp-125L,
    -0x8671cb6dbfc294a3p-133L,
};
inline constexpr ld k_sinh[] = {
    0xaaaaaaaaaaaaaaabp-66L,
    0x8888888888888889p-70L,
    0xd00d00d00d00d00dp-76L,
    0xb8ef1d2ab6399c7dp-82L,
    0xd7322b3faa271c7fp-89L,
    0xb092309d43684be5p-96L,
    0xd73f9f399dc0f88fp-104L,
    0xca963b81856a5359p-112L,
    0x97a4da340a0ab926p-120L,
    0xb8dc77b6e7ab8c5fp-129L,
};
inline constexpr ld k_sinh_tail[] = {
    0x8888888888888889p-70L,
    0xd00d00d00d00d00dp-76L,
    0xb8ef1d2ab6399c7dp-82L,
    0xd7322b3faa271c7fp-89L,
    0xb092309d43684be5p-96L,
    0xd73f9f399dc0f88fp-104L,
    0xca963b81856a5359p-112L,
    0x97a4da340a0ab926p-120L,
    0xb8dc77b6e7ab8c5fp-129L,
};
inline constexpr ld k_exp[] = {
    0x8000000000000000p-64L,
    0xaaaaaaaaaaaaaaabp-66L,
    0xaaaaaaaaaaaaaaabp-68L,
    0x8888888888888889p-70L,
    0xb60b60b60b60b60bp-73L,
    0xd00d00d00d00d00dp-76L,
    0xd00d00d00d00d00dp-79L,
    0xb8ef1d2ab6399c7dp-82L,
    0x93f27dbbc4fae397p-85L,
    0xd7322b3faa271c7fp-89L,
    0x8f76c77fc6c4bdaap-92L,
    0xb092309d43684be5p-96L,
    0xc9cba54603e4e906p-100L,
    0xd73f9f399dc0f88fp-104L,
    0xd73f9f399dc0f88fp-108L,
    0xca963b81856a5359p-112L,
    0xb413c31dcbecbbdep-116L,
};
inline constexpr ld k_exp_tail[] = {
    0xaaaaaaaaaaaaaaabp-66L,
    0xaaaaaaaaaaaaaaabp-68L,
    0x8888888888888889p-70L,
    0xb60b60b60b60b60bp-73L,
    0xd00d00d00d00d00dp-76L,
    0xd00d00d00d00d00dp-79L,
    0xb8ef1d2ab6399c7dp-82L,
    0x93f27dbbc4fae397p-85L,
    0xd7322b3faa271c7fp-89L,
    0x8f76c77fc6c4bdaap-92L,
    0xb092309d43684be5p-96L,
    0xc9cba54603e4e906p-100L,
    0xd73f9f399dc0f88fp-104L,
    0xd73f9f399dc0f88fp-108L,
    0xca963b81856a5359p-112L,
    0xb413c31dcbecbbdep-116L,
};
inline constexpr ld k_log_tail[] = {
    0xe38e38e38e38e38ep-67L,
    0xba2e8ba2e8ba2e8cp-67L,
    0x9d89d89d89d89d8ap-67L,
    0x8888888888888889p-67L,
    0xf0f0f0f0f0f0f0f1p-68L,
    0xd79435e50d79435ep-68L,
    0xc30c30c30c30c30cp-68L,
    0xb21642c8590b2164p-68L,
    0xa3d70a3d70a3d70ap-68L,
    0x97b425ed097b425fp-68L,
    0x8d3dcb08d3dcb08dp-68L,
    0x8421084210842108p-68L,
    0xf83e0f83e0f83e10p-69L,
    0xea0ea0ea0ea0ea0fp-69L,
};
inline constexpr ld k_lg1[] = {
    -0xcd26aadf5596b7adp-65L,
    0x8a8991563ec241b6p-65L,
    -0xd45ce0bd530a492fp-66L,
    0xada06588061830a5p-66L,
    -0x9381d0ee751d72c6p-66L,
    0x80859b57c31cb746p-66L,
    -0xe4033836abeda2adp-67L,
    0xcd00f1c2eafc7981p-67L,
    -0xba461988a636b07bp-67L,
    0xaab56b1921de202cp-67L,
    -0x9d8ecb8fe2cc261ap-67L,
    0x924b6fc19062fb9fp-67L,
    -0x88899a3b73ff01ddp-67L,
    0x8000803266f59178p-67L,
    -0xf0f16988f455e4e6p-68L,
    0xe38e71d105a0c483p-68L,
    -0xd79450dab44502e9p-68L,
    0xccccd99a96add052p-68L,
    -0xc30c36dbbdfdfadcp-68L,
    0xba2e8e8bbc6fc881p-68L,
    -0xb216442c8db365dfp-68L,
    0xaaaaab5557ee6aa9p-68L,
    -0xa3d70a8f5cfdbba0p-68L,
    0x9d89d8c4ec92f3e9p-68L,
    -0x97b426000015e20bp-68L,
    0x9249249b6dbde3e5p-68L,
    -0x8d3dcb0d3dcd4c3ep-68L,
    0x8888888aaaab655bp-68L,
    -0x8421084318c66dc5p-68L,
    0x8000000080001372p-68L,
    -0xf83e0f845d175263p-69L,
};
inline constexpr ld k_lg2[] = {
    -0x89f000d2abb03409p-67L,
    0xa8991563ec241b60p-69L,
    -0xf2027e10c7af8c37p-71L,
    0xbd6eb756db617ea5p-72L,
    -0x9c562e15fc703e76p-73L,
    0x859b57c31cb745f3p-74L,
    -0xe9fea63b697e3e38p-76L,
    0xd093d878beb2d19dp-77L,
    -0xbc6f2debe40f7797p-78L,
    0xac06e77337581126p-79L,
    -0x9e5e4b1e7112142bp-80L,
    0x92cbd1cf9a555c81p-81L,
    -0x88d975bb3caa08e4p-82L,
    0x803266f5917879d0p-83L,
    -0xf13006c9e7e975dap-85L,
    0xe3b5dd9f83d26bb3p-86L,
    -0xd7ad365dfc54bb2cp-87L,
    0xccdc9e1038587a06p-88L,
    -0xc31639a6f9f56366p-89L,
    0xba34ed667d6e6593p-90L,
    -0xb21a54223d75681bp-91L,
    0xaaad43bffe9614f1p-92L,
    -0xa3d8b3c92c68720ap-93L,
};
inline constexpr ld k_stirling[] = {
    0xaaaaaaaaaaaaaaabp-67L,
    -0xb60b60b60b60b60bp-72L,
    0xd00d00d00d00d00dp-74L,
    -0x9c09c09c09c09c0ap-74L,
    0xdca8f158c7f91ab8p-74L,
    -0xfb5586ccc9e3e410p-73L,
    0xd20d20d20d20d20dp-71L,
    -0xf21436587a9cbee1p-69L,
    0xb7f4b1c0f033ffd1p-66L,
    -0xb23b3808c0f9cf6ep-63L,
    0xd672219167002d3ap-60L,
};
inline constexpr ld k_stirling_tail[] = {
    -0xb60b60b60b60b60bp-72L,
    0xd00d00d00d00d00dp-74L,
    -0x9c09c09c09c09c0ap-74L,
    0xdca8f158c7f91ab8p-74L,
    -0xfb5586ccc9e3e410p-73L,
    0xd20d20d20d20d20dp-71L,
    -0xf21436587a9cbee1p-69L,
    0xb7f4b1c0f033ffd1p-66L,
    -0xb23b3808c0f9cf6ep-63L,
    0xd672219167002d3ap-60L,
};
inline constexpr std::uint64_t k_two_over_pi[262] = {
    0xa2f9836e4e441529u, 0xfc2757d1f534ddc0u, 0xdb6295993c439041u, 0xfe5163abdebbc561u,
    0xb7246e3a424dd2e0u, 0x06492eea09d1921cu, 0xfe1deb1cb129a73eu, 0xe88235f52ebb4484u,
    0xe99c7026b45f7e41u, 0x3991d639835339f4u, 0x9c845f8bbdf9283bu, 0x1ff897ffde05980fu,
    0xef2f118b5a0a6d1fu, 0x6d367ecf27cb09b7u, 0x4f463f669e5fea2du, 0x7527bac7ebe5f17bu,
    0x3d0739f78a5292eau, 0x6bfb5fb11f8d5d08u, 0x56033046fc7b6babu, 0xf0cfbc209af4361du,
    0xa9e391615ee61b08u, 0x6599855f14a06840u, 0x8dffd8804d732731u, 0x06061556ca73a8c9u,
    0x60e27bc08c6b47c4u, 0x19c367cddce8092au, 0x8359c4768b961ca6u, 0xddaf44d15719053eu,
    0xa5ff07053f7e33e8u, 0x32c2de4f98327dbbu, 0xc33d26ef6b1e5ef8u, 0x9f3a1f35caf27f1du,
    0x87f121907c7c246au, 0xfa6ed5772d30433bu, 0x15c614b59d19c3c2u, 0xc4ad414d2c5d000cu,
    0x467d862d71e39ac6u, 0x9b0062337cd2b497u, 0xa7b4d55537f63ed7u, 0x1810a3fc764d2a9du,
    0x64abd770f87c6357u, 0xb07ae715175649c0u, 0xd9d63b3884a7cb23u, 0x24778ad623545ab9u,
    0x1f001b0af1dfce19u, 0xff319f6a1e666157u, 0x9947fbacd87f7eb7u, 0x652289e83260bfe6u,
    0xcdc4ef09366cd43fu, 0x5dd7de16de3b5892u, 0x9bde2822d2e88628u, 0x4d58e232cac616e3u,
    0x08cb7de050c017a7u, 0x1df35be01834132eu, 0x6212830148835b8eu, 0xf57fb0adf2e91e43u,
    0x4a48d36710d8ddaau, 0x425faece616aa428u, 0x0ab499d3f2a6067fu, 0x775c83c2a3883c61u,
    0x78738a5a8cafbdd7u, 0x6f63a62dcbbff4efu, 0x818d67c12645ca55u, 0x36d9cad2a8288d61u,
    0xc277c9121426049bu, 0x4612c459c444c5c8u, 0x91b24df31700ad43u, 0xd4e5492910d5fdfcu,
    0xbe00cc941eeece70u, 0xf53e1380f1ecc3e7u, 0xb328f8c79405933eu, 0x71c1b3092ef3450bu,
    0x9c12887b20ab9fb5u, 0x2ec292472f327b6du, 0x550c90a7721fe76bu, 0x96cb314a1679e279u,
    0x4189dff49794e884u, 0xe6e29731996bed88u, 0x365f5f0efdbbb49au, 0x486ca46742727132u,
    0x5d8db8159f09e5bcu, 0x25318d3974f71c05u, 0x30010c0d68084b58u, 0xee2c90aa4702e774u,
    0x24d6bda67df77248u, 0x6eef169fa6948ef6u, 0x91b45153d1f20acfu, 0x3398207e4bf56863u,
    0xb25f3edd035d407fu, 0x8985295255c06437u, 0x10d86d324832754cu, 0x5bd4714e6e5445c1u,
    0x090b69f52ad56614u, 0x9d072750045ddb3bu, 0xb4c576ea17f9877du, 0x6b49ba271d296996u,
    0xacccc65414ad6ae2u, 0x9089d98850722cbeu, 0xa4049407777030f3u, 0x27fc00a871ea49c2u,
    0x663de06483dd9797u, 0x3fa3fd94438c860du, 0xde41319d39928c70u, 0xdde7b7173bdf082bu,
    0x3715a0805c93805au, 0x921110d8e80faf80u, 0x6c4bffdb0f903876u, 0x185915a562bbcb61u,
    0xb989c7bd401004f2u, 0xd2277549f6b6ebbbu, 0x22dbaa140a2f2689u, 0x768364333b091a94u,
    0x0eaa3a51c2a31daeu, 0xedaf12265c4dc26du, 0x9c7a2d9756c0833fu, 0x03f6f0098c402b99u,
    0x316d07b43915200cu, 0x5bc3d8c492f54badu, 0xc6a5ca4ecd37a736u, 0xa9e69492ab6842ddu,
    0xde6319ef8c76528bu, 0x6837dbfcaba1ae31u, 0x15dfa1ae00dafb0cu, 0x664d64b705ed3065u,
    0x29bf56573aff47b9u, 0xf96af3be75df9328u, 0x3080abf68c6615cbu, 0x040622fa1de4d9a4u,
    0xb33d8f1b5709cd36u, 0xe9424ea4be13b523u, 0x331aaaf0a8654fa5u, 0xc1d20f3f0bcd785bu,
    0x76f923048b7b7217u, 0x8953a6c6e26e6f00u, 0xebef584a9bb7dac4u, 0xba66aacfcf761d02u,
    0xd12df1b1c1998c77u, 0xadc3da4886a05df7u, 0xf480c62ff0ac9aecu, 0xddbc5c3f6dded01fu,
    0xc790b6db2a3a25a3u, 0x9aaf009353ad0457u, 0xb6b42d297e804ba7u, 0x07da0eaa76a1597bu,
    0x2a12162db7dcfde5u, 0xfafedb89fdbe896cu, 0x76e4fca90670803eu, 0x156e85ff87fd073eu,
    0x2833676186182aeau, 0xbd4dafe7b36e6d8fu, 0x3967955bbf3148d7u, 0x8416df30432dc735u,
    0x6125ce70c9b8cb30u, 0xfd6cbfa200a4e46cu, 0x05a0dd5a476f21d2u, 0x1262845cb9496170u,
    0xe0566b0152993755u, 0x50b7d51ec4f1335fu, 0x6e13e4305da92e85u, 0xc3b21d3632a1a4b7u,
    0x08d4b1ea21f716e4u, 0x698f77ff2780030cu, 0x2d408da0cd4f99a5u, 0x20d3a2b30a5d2f42u,
    0xf9b4cbda11d0be7du, 0xc1db9bbd17ab81a2u, 0xca5c6a0817552e55u, 0x0027f0147f8607e1u,
    0x640b148d4196debeu, 0x872afddab6256b34u, 0x897bfef3059ebfb9u, 0x4f6a68a82a4a5ac4u,
    0x4fbcf82d985ad795u, 0xc7f48d4d0da63a20u, 0x5f57a4b13f149538u, 0x800120cc86dd71b6u,
    0xdec9f560bf11654du, 0x6b0701acb08cd0c0u, 0xb24855510efb1ec3u, 0x72953b06a33540c0u,
    0x7bdc06cc45e0fa29u, 0x4ec8cad641f3e8deu, 0x647cd8649b31bed9u, 0xc397a4d45877c5e3u,
    0x6913daf03c3aba46u, 0x18465f7555f5bdd2u, 0xc6926e5d2eaced44u, 0x0e423e1c87c461e9u,
    0xfd29f3d6e7ca7c22u, 0x35916fc5e0088dd7u, 0xffe26a6ec6fdb0c1u, 0x0893745d7cb2ad6bu,
    0x9d6ecd7b723e6a11u, 0xc6a9cff7df7329bau, 0xc9b55100b70db2e2u, 0x24ba74607de58ad8u,
    0x742c150d0c188194u, 0x667e162901767a9fu, 0xbefdfdef4556367eu, 0xd913d9ecb9ba8bfcu,
    0x97c427a831c36ef1u, 0x36c59456a8d8b5a8u, 0xb40ecccf2d891234u, 0x576f89562ce3ce99u,
    0xb920d6aa5e6b9c2au, 0x3ecc5f114a0bfdfbu, 0xf4e16d3b8e2c86e2u, 0x84d4e9a9b4fcd1eeu,
    0xefc9352e61392f44u, 0x2138c8d91b0afc81u, 0x6a4afbd81c2f84b4u, 0x538c994ecc2254dcu,
    0x552ad6c6c096190bu, 0xb8701a649569605au, 0x26ee523f0f117f11u, 0xb5f4f5cbfc2dbc34u,
    0xeebc34cc5de8605eu, 0xdd9b8e67ef3392b8u, 0x17c99b5861bc57e1u, 0xc68351103ed84871u,
    0xdddd1c2da118af46u, 0x2c21d7f359987ad9u, 0xc0549efa864ffc06u, 0x56ae79e536228922u,
    0xad38dc9367aae855u, 0x3826829be7caa40du, 0x51b133990ed7a948u, 0x0569f0b265a7887fu,
    0x974c8836d1f9b392u, 0x214a827b21cf98dcu, 0x9f405547dc3a74e1u, 0x42eb67df9dfe5fd4u,
    0x5ea4677b7aacbaa2u, 0xf65523882b55ba41u, 0x086e59862a218347u, 0x39e6e389d49ee540u,
    0xfb49e956ffca0f1cu, 0x8a59c52bfa94c5c1u, 0xd3cfc50fae5adb86u, 0xc5476243853b8621u,
    0x94792c8761107b4cu, 0x2a1a2c8012bf4390u, 0x2688893c78e4c4a8u, 0x7bdbe5c23ac4eaf4u,
    0x268a67f7bf920d2bu, 0xa365b1933d0b7cbdu, 0xdc51a463dd27dde1u, 0x6919949a9529a828u,
    0xce68b4ed09209f44u, 0xca984e638270237cu, 0x7e32b90f8ef5a7e7u, 0x561408f1212a9db5u,
    0x4d7e6f5119a5abf9u, 0xb5d6df8261dd9602u, 0x36169f3ac4a1a283u, 0x6ded727a8d39a9b8u,
    0x825c326b5b2746edu, 0x34007700d255f4fcu, 0x4d59018071e0e13fu, 0x89b295f364a8f1aeu,
    0xa74b38fc4ceab2bbu, 0x47270babc3a734bau,
};

// ---- trigonometric functions ----------------------------------------------------
// q, r with |x| = q pi/2 + r (q mod 4), |r| <= pi/4 (to rounding), r double-length.
// |x| >= pi/4, finite.  x = m 2^e with m a 64-bit integer; the bits of 2/pi worth less than
// 4 / (m 2^e) contribute multiples of 4 and are skipped; 256 further bits give the fraction of
// x 2/pi to within 2^(64 + 1 - 255) = 2^-190 (absolute, in quadrants).
[[nodiscard]] inline int rem_pio2(ld ax, dd& r) noexcept {
    int ex;
    const ld f = std::frexp(ax, &ex);                                  // ax = f 2^ex, f in [0.5, 1)
    const std::uint64_t m = static_cast<std::uint64_t>(std::ldexp(f, 64));  // exact
    const int e = ex - 64;
    const int i0 = e - 1 > 1 ? e - 1 : 1;                              // first bit of 2/pi used
    std::uint64_t w[4];
    const int k0 = (i0 - 1) / 64, sh = (i0 - 1) % 64;
    for (int j = 0; j < 4; ++j)
        w[j] = sh == 0 ? k_two_over_pi[k0 + j]
                       : (k_two_over_pi[k0 + j] << sh) | (k_two_over_pi[k0 + j + 1] >> (64 - sh));
    // P = m * W, 320 bits in p[0] (most significant) .. p[4]
    std::uint64_t p[5];
    u128 carry = 0;
    for (int j = 3; j >= 0; --j) {
        const u128 t = static_cast<u128>(m) * w[j] + carry;
        p[j + 1] = static_cast<std::uint64_t>(t);
        carry = t >> 64;
    }
    p[0] = static_cast<std::uint64_t>(carry);
    // value(P) = m W 2^(e - i0 - 255): pt fraction bits, pt in [254, 320]
    const int pt = i0 + 255 - e;
    auto bit = [&](int i) -> unsigned {  // bit i of P, i = 0 the least significant
        return static_cast<unsigned>((p[4 - i / 64] >> (i % 64)) & 1u);
    };
    int q = static_cast<int>((pt < 320 ? bit(pt) : 0u) | ((pt + 1 < 320 ? bit(pt + 1) : 0u) << 1));
    // F = P mod 2^pt: clear the integer bits
    for (int i = 0; i < 5; ++i) {
        const int lo_bit = (4 - i) * 64;  // index of the least significant bit of p[i]
        if (lo_bit >= pt) p[i] = 0;
        else if (lo_bit + 64 > pt) p[i] &= (~std::uint64_t(0)) >> (64 - (pt - lo_bit));
    }
    bool negative = false;
    if (bit(pt - 1)) {  // fraction >= 1/2: next quadrant, F = 2^pt - F
        q += 1;
        negative = true;
        std::uint64_t c = 1;
        for (int i = 4; i >= 0; --i) {
            const std::uint64_t v = ~p[i] + c;
            c = (c && v == 0) ? 1 : 0;
            p[i] = v;
        }
        for (int i = 0; i < 5; ++i) {  // keep pt bits
            const int lo_bit = (4 - i) * 64;
            if (lo_bit >= pt) p[i] = 0;
            else if (lo_bit + 64 > pt) p[i] &= (~std::uint64_t(0)) >> (64 - (pt - lo_bit));
        }
    }
    int L = -1;  // leading one
    for (int i = 0; i < 5 && L < 0; ++i)
        if (p[i]) L = (4 - i) * 64 + 63 - __builtin_clzll(p[i]);
    if (L < 0) {
        r = {0.0L, 0.0L};
        return q & 3;
    }
    auto word = [&](int i) -> std::uint64_t { return i >= 0 && i < 5 ? p[4 - i] : 0; };  // bits 64i ..
    auto bits64 = [&](int top) -> std::uint64_t {  // bits top-63 .. top of P (zeros below 0)
        const int low = top - 63;
        if (low < 0) return word(0) << (-low);  // top < 63
        const int i = low / 64, sh = low % 64;
        return sh == 0 ? word(i) : (word(i) >> sh) | (word(i + 1) << (64 - sh));
    };
    const std::uint64_t h = bits64(L), l = bits64(L - 64);
    const dd fr = fast_two_sum(static_cast<ld>(h), std::ldexp(static_cast<ld>(l), -64));
    dd rr = mul(scale(fr, L - 63 - pt), k_pio2);
    if (negative) rr = neg(rr);
    r = rr;
    return q & 3;
}
// Cody-Waite for pi/4 < |x| < 2^23 pi/2 (k < 2^23): r = x - k (c1 + c2 + c3 + c4), pi/2 to
// 232 bits with c1 of 40 bits, so x - k c1 is exact and r has an absolute error below 2^-200
// (the smallest |r| of a long double in this range is far larger).
[[nodiscard]] inline int rem_pio2_cw(ld ax, dd& r) noexcept {
    const ld kf = round_int(ax * k_2_over_pi);
    const ld t = ax - kf * k_pio2_cw[0];  // exact
    const dd p2 = two_prod(kf, k_pio2_cw[1]), p3 = two_prod(kf, k_pio2_cw[2]);
    dd s = two_sum(t, -p2.hi);
    s = add(s, -p2.lo);
    s = add(s, neg(p3));
    r = add(s, -kf * k_pio2_cw[3]);
    return static_cast<int>(static_cast<long long>(kf) & 3);
}
// sin and cos of r = hi + lo, |r| <= pi/4 (+ rounding), double-length results.
[[nodiscard]] inline dd ksin(dd r) noexcept {
    const dd z = two_prod(r.hi, r.hi);
    const ld P = k_sin[0] + z.hi * horner(k_sin_tail, z.hi);  // -1/6 + r^2/120 - ...
    const dd t = mul(mul(z, P), r.hi);                         // r^3 P
    return add(fast_two_sum(r.hi, t.hi), t.lo + r.lo * (1.0L - 0.5L * z.hi));
}
[[nodiscard]] inline dd kcos(dd r) noexcept {
    const dd z = two_prod(r.hi, r.hi);
    const ld hz = 0.5L * z.hi, w = 1.0L - hz;
    const ld tail = (((1.0L - w) - hz) - 0.5L * z.lo) + (z.hi * z.hi * horner(k_cos, z.hi) - r.hi * r.lo);
    return fast_two_sum(w, tail);
}
// sin r or cos r (|r| <= pi/4) to full double-length accuracy (about 2^-120), by the Taylor
// series in dd arithmetic.  ksin/kcos above are accurate to about 2^-70, enough for a rounded
// result; the gamma reflection needs log|sin(pi x)| next to the zeros of lgamma.
[[nodiscard]] inline dd trig_series(dd r, bool cosine) noexcept {
    const dd z = neg(mul(r, r));
    dd term = cosine ? dd{1.0L, 0.0L} : r, sum = term;
    for (int k = cosine ? 1 : 2; k < 40; k += 2) {  // term *= -r^2 / (k (k + 1))
        term = div(mul(term, z), {static_cast<ld>(k * (k + 1)), 0.0L});
        sum = add(sum, term);
        if (std::fabs(term.hi) < std::fabs(sum.hi) * 0x1p-130L) break;
    }
    return sum;
}
// sin(pi a) for |a| <= 1/2, double-length (for the gamma reflection)
[[nodiscard]] inline dd sinpi_dd(ld a) noexcept {
    const ld b = std::fabs(a);
    const dd s = b <= 0.25L ? trig_series(mul(k_pi, b), false) : trig_series(mul(k_pi, 0.5L - b), true);
    return a < 0 ? neg(s) : s;
}
enum class trig { sin, cos, tan };
[[nodiscard]] inline ld trig_eval(ld x, trig f) noexcept {
    const ld ax = std::fabs(x);
    if (ax < 0x1p-32L) return f == trig::cos ? 1.0L : x;  // also +-0, and subnormals
    dd r{ax, 0.0L};
    int q = 0;
    if (ax > k_pio4) q = ax < 0x1.92p+23L ? rem_pio2_cw(ax, r) : rem_pio2(ax, r);
    ld res;
    if (f == trig::tan) {
        const dd s = ksin(r), c = kcos(r);
        res = (q & 1) ? -round(div(c, s)) : round(div(s, c));
    } else {
        const int qq = f == trig::cos ? (q + 1) & 3 : q;  // cos x = sin(x + pi/2)
        const dd v = (qq & 1) ? kcos(r) : ksin(r);
        res = (qq & 2) ? -round(v) : round(v);
    }
    return (f != trig::cos && x < 0) ? -res : res;
}

// ---- exp, expm1, sinh, log, pow -------------------------------------------------
// e^r - 1 for r = hi + lo, |r| <= ln2/2 (+ rounding), double-length.
[[nodiscard]] inline dd kexpm1(dd r) noexcept {
    const dd z = two_prod(r.hi, r.hi);
    const ld P = k_exp[0] + r.hi * horner(k_exp_tail, r.hi);  // 1/2 + r/6 + ...
    const dd t = mul(z, P);
    return add(fast_two_sum(r.hi, t.hi), t.lo + r.lo * (1.0L + r.hi));
}
// y = k ln2 + r, |r| <= ln2/2 (+ rounding), for |y.hi| < 2^15 ln2 (about 22713).
[[nodiscard]] inline int exp_reduce(dd y, dd& r) noexcept {
    const ld kf = round_int(y.hi * k_inv_ln2);
    const dd a = add(two_sum(y.hi - kf * k_ln2_1, y.lo), neg(two_prod(kf, k_ln2_2)));
    r = add(a, -kf * k_ln2_3);
    return static_cast<int>(kf);
}
// e^y rounded; y double-length
[[nodiscard]] inline ld exp_dd(dd y) noexcept {
    if (y.hi > 11400.0L) return std::numeric_limits<ld>::infinity();
    if (y.hi < -11450.0L) return 0.0L;
    dd r;
    const int k = exp_reduce(y, r);
    const dd v = add(kexpm1(r), 1.0L);  // e^r
    return scale(round(v), k);
}
[[nodiscard]] inline ld expm1(ld x) noexcept {
    if (!(std::fabs(x) <= 11400.0L) || x == 0) return std::expm1(x);  // NaN, +-inf, +-0, overflow
    if (std::fabs(x) < 0x1p-65L) return x;
    if (x < -46.0L) return -1.0L;  // e^x < 2^-66
    if (std::fabs(x) <= 0.34657359027997265471L) return round(kexpm1({x, 0.0L}));
    dd r;
    const int k = exp_reduce({x, 0.0L}, r);
    const dd em = kexpm1(r);
    if (k > 64) return scale(round(add(em, 1.0L)), k);  // the -1 is below the last bit
    return round(add(scale(em, k), two_sum(pow2(k), -1.0L)));  // 2^k em + (2^k - 1)
}
[[nodiscard]] inline ld sinh(ld x) noexcept {
    const ld ax = std::fabs(x);
    if (!(ax <= 11400.0L) || x == 0) return std::sinh(x);
    if (ax < 0x1p-32L) return x;
    ld res;
    if (ax <= 1.0L) {  // x + x^3 (1/6 + x^2/120 + ...)
        const dd z = two_prod(ax, ax);
        const ld P = k_sinh[0] + z.hi * horner(k_sinh_tail, z.hi);
        const dd t = mul(mul(z, P), ax);
        res = round(add(fast_two_sum(ax, t.hi), t.lo));
    } else {
        dd r;
        const int k = exp_reduce({ax, 0.0L}, r);
        const dd v = add(kexpm1(r), 1.0L);  // e^ax = 2^k v
        if (ax > 46.0L) res = scale(round(v), k - 1);
        else res = 0.5L * round(add(scale(v, k), neg(scale(div({1.0L, 0.0L}, v), -k))));
    }
    return x < 0 ? -res : res;
}
// log x, double-length, x > 0 finite: x = 2^e m, m in [sqrt(1/2), sqrt(2)),
// log m = 2 atanh(f), f = (m - 1)/(m + 1), |f| <= 0.1716, s = f^2 <= 0.0295.
[[nodiscard]] inline dd log_dd(ld x) noexcept {
    int e;
    ld m;
    if (const int be = biased_exponent(x); be != 0) {  // normal: x = m 2^e, m in [0.5, 1)
        e = be - 16382;
        m = scale(x, -e);  // exact
    } else {
        m = std::frexp(x, &e);
    }
    if (m < 0.70710678118654752440L) m *= 2, --e;
    const dd f = div({m - 1.0L, 0.0L}, two_sum(m, 1.0L));
    const dd s = mul(f, f);
    // 1 + s/3 + s^2/5 + s^3/7 double-length, s^4 (1/9 + s/11 + ...) in long double (< 2^-23)
    const dd s2 = mul(s, s), s3 = mul(s2, s);
    dd P = add(add(mul(s, k_inv3), mul(s2, k_inv5)), mul(s3, k_inv7));
    P = add(add(P, s2.hi * s2.hi * horner(k_log_tail, s.hi)), 1.0L);
    const dd lm = mul(mul(f, P), 2.0L);
    const ld ef = static_cast<ld>(e);
    const dd le = add(add(two_prod(ef, k_ln2_2), ef * k_ln2_3), ef * k_ln2_1);  // ef*k_ln2_1 exact
    return add(le, lm);
}
[[nodiscard]] inline dd log_dd(dd x) noexcept { return add(log_dd(x.hi), x.lo / x.hi); }
[[nodiscard]] inline ld pow_pos(ld x, ld y) noexcept {  // x > 0, x != 1, both finite
    const dd l = log_dd(x);
    const ld est = l.hi * y;  // overflow and underflow decided before the exact product
    if (est > 11400.0L) return std::numeric_limits<ld>::infinity();
    if (est < -11450.0L) return 0.0L;
    const dd t = mul(l, y);
    return exp_dd(t);
}
[[nodiscard]] inline ld pow(ld x, ld y) noexcept {
    if (!std::isfinite(x) || !std::isfinite(y) || x == 0 || y == 0 || std::fabs(x) == 1.0L)
        return std::pow(x, y);
    if (x > 0) return pow_pos(x, y);
    if (std::nearbyint(y) != y) return std::pow(x, y);  // NaN (negative base, non-integer y)
    const ld r = pow_pos(-x, y);
    return std::fmod(y, 2.0L) != 0 ? -r : r;
}

// ---- lgamma, tgamma ---------------------------------------------------------------
// psi(y) to about 1e-7 relative (only multiplies the tiny rounding error of a shifted argument)
[[nodiscard]] inline ld digamma_approx(ld y) noexcept {
    ld acc = 0;
    for (; y < 10.0L; y += 1.0L) acc -= 1.0L / y;
    const ld w2 = 1.0L / (y * y);
    return acc + std::log(y) - 0.5L / y - w2 * (1.0L / 12 - w2 * (1.0L / 120 - w2 / 252));
}
// log|Gamma(x)| for x >= 20: Stirling with 11 terms (error < 2^-80 there).
[[nodiscard]] inline dd lgamma_stirling(ld x) noexcept {
    // 1/(12 x) double-length, then w^3 (B4/12 + B6/30 w^2 + ...) in long double
    const ld w = 1.0L / x, w2 = w * w;
    const dd series = add(div({1.0L, 0.0L}, two_prod(12.0L, x)), w * w2 * horner(k_stirling_tail, w2));
    dd t = mul(two_sum(x, -0.5L), log_dd(x));
    t = add(add(t, -x), k_hl2pi);
    return add(t, series);
}
// log Gamma(x) for x > 0 finite, double-length (absolute error about 2^-80; relative next to
// the zeros at 1 and 2, where the Taylor series are used).
[[nodiscard]] inline dd lgamma_pos(ld x) noexcept {
    if (std::fabs(x - 1.0L) <= 0.2L) {  // lgamma(1 + z) = -gamma z + zeta(2)/2 z^2 - ...
        const ld z = x - 1.0L;             // exact
        dd s = add(mul(k_egamma, -z), mul(two_prod(z, z), k_zeta2h.hi));
        s = add(s, z * z * k_zeta2h.lo + z * z * z * horner(k_lg1, z));
        return s;
    }
    if (std::fabs(x - 2.0L) <= 0.2L) {  // lgamma(2 + z) = (1 - gamma) z + (zeta(2) - 1)/2 z^2 - ...
        const ld z = x - 2.0L;
        dd s = add(mul(k_1mgamma, z), mul(two_prod(z, z), k_zeta2m1h.hi));
        s = add(s, z * z * k_zeta2m1h.lo + z * z * z * horner(k_lg2, z));
        return s;
    }
    if (x > 0x1p62L) {  // x (log x - 1) + (1/2 log(2 pi) - 1/2 log x), with x = m 2^e kept unscaled
        const int e = biased_exponent(x) - 16383;
        const ld m = scale(x, -e);  // exact, in [1, 2)
        const dd L = log_dd(x);
        dd s = mul(add(L, -1.0L), m);
        if (e < 16000) s = add(s, scale(add(mul(L, -0.5L), k_hl2pi), -e));
        return {scale(round(s), e), 0.0L};
    }
    if (x >= 20.0L) return lgamma_stirling(x);
    // lgamma(x) = lgamma(x + n) - log(x (x + 1) ... (x + n - 1)), x + n >= 20
    const int n = static_cast<int>(std::ceil(20.0L - x));
    dd prod{x, 0.0L};
    for (int i = 1; i < n; ++i) prod = mul(prod, two_sum(x, static_cast<ld>(i)));
    const dd y = two_sum(x, static_cast<ld>(n));  // x + n, exactly
    dd lg = lgamma_stirling(y.hi);
    // y.lo psi(y.hi), psi(y) ~ log y - 1/(2y) (y >= 20, |y.lo| <= 2^-59)
    lg = add(lg, y.lo * digamma_approx(y.hi));
    return add(lg, neg(log_dd(prod)));
}
// x < 0, not an integer: lgamma(x) = log pi - log|sin(pi x)| - lgamma(1 - x); sign of Gamma.
[[nodiscard]] inline dd lgamma_neg(ld x, int& sign) noexcept {
    const ld n = std::nearbyint(x), a = x - n;  // exact, |a| <= 1/2
    const dd s = sinpi_dd(a);                   // sin(pi x) = (-1)^n sin(pi a)
    const bool odd = std::fmod(n, 2.0L) != 0;
    sign = ((s.hi < 0) != odd) ? -1 : 1;
    const dd y = two_sum(1.0L, -x);  // 1 - x
    dd lg = lgamma_pos(y.hi);
    if (y.lo != 0) lg = add(lg, y.lo * digamma_approx(y.hi));
    const dd ls = log_dd(s.hi < 0 ? neg(s) : s);
    return add(add(k_logpi, neg(ls)), neg(lg));
}
[[nodiscard]] inline ld lgamma(ld x) noexcept {
    if (!std::isfinite(x) || x == 0 || (x < 0 && std::nearbyint(x) == x)) return std::lgamma(x);  // poles
    if (std::fabs(x) < 0x1p-64L) return -round(log_dd(std::fabs(x)));  // - log|x| - gamma x
    if (x > 0) return round(lgamma_pos(x));
    int sign;
    return round(lgamma_neg(x, sign));
}
[[nodiscard]] inline ld tgamma(ld x) noexcept {
    if (!std::isfinite(x) || x == 0 || (x < 0 && std::nearbyint(x) == x) || x > 1756.0L || x < -0x1p62L)
        return std::tgamma(x);
    if (std::fabs(x) < 0x1p-16380L) return 1.0L / x;  // gamma x is below 2^-16380 of it
    if (std::fabs(x) < 0x1p-64L) return round(add(div({1.0L, 0.0L}, {x, 0.0L}), -k_egamma.hi));  // 1/x - gamma
    if (x > 0) {
        if (x == std::nearbyint(x) && x <= 30.0L) {  // (x - 1)!, exact while it fits
            ld f = 1.0L;
            bool exact = true;
            for (int i = 2; i < static_cast<int>(x) && exact; ++i) {
                const dd p = two_prod(f, static_cast<ld>(i));
                exact = p.lo == 0;
                f = p.hi;
            }
            if (exact) return f;
        }
        return exp_dd(lgamma_pos(x));
    }
    int sign;
    const ld g = exp_dd(lgamma_neg(x, sign));
    return sign < 0 ? -g : g;
}
// sin x and cos x (x finite) with one reduction, each as lm::sin/lm::cos would give
inline void sincos(ld x, ld& s, ld& c) noexcept {
    const ld ax = std::fabs(x);
    if (ax < 0x1p-32L) {
        s = x;
        c = 1.0L;
        return;
    }
    dd r{ax, 0.0L};
    int q = 0;
    if (ax > k_pio4) q = ax < 0x1.92p+23L ? rem_pio2_cw(ax, r) : rem_pio2(ax, r);
    const ld sr = round(ksin(r)), cr = round(kcos(r));
    switch (q) {
        case 0: s = sr, c = cr; break;
        case 1: s = cr, c = -sr; break;
        case 2: s = -sr, c = -cr; break;
        default: s = -cr, c = sr; break;
    }
    if (x < 0) s = -s;
}
// asin and acos for 1/2 < |x| <= 1, where Apple's asinl and acosl lose up to about 1100 ulp
// next to +-1 (fuzz, 1 - |x| down to 6e-12).  t = sqrt((1 - a)(1 + a)), a = |x|, double-length
// (1 - a exact, 1 + a and the product as dd, then one Newton step for the square root); then
// asin a = atan2(a, t) and acos x = atan2(t, x), with the first-order correction for t.lo
// (d/dt atan2(a, t) = -a and d/dt atan2(t, x) = x, as a^2 + t^2 = 1).  atan2l is accurate on
// Apple's libm (the same bits as glibc's) and is well conditioned there.
inline void asin_acos_t(ld a, ld& th, ld& tl) noexcept {
    const dd p = mul(fast_two_sum(1.0L, a), 1.0L - a);
    th = std::sqrt(p.hi);
    const dd sq = two_prod(th, th);
    tl = th == 0 ? 0.0L : ((p.hi - sq.hi) - sq.lo + p.lo) / (2 * th);
}
[[nodiscard]] inline ld asin(ld x) noexcept {
    const ld a = std::fabs(x);
    if (!(a > 0.5L) || a > 1) return std::asin(x);
    ld th, tl;
    asin_acos_t(a, th, tl);
    return std::copysign(std::atan2(a, th) - a * tl, x);
}
[[nodiscard]] inline ld acos(ld x) noexcept {
    const ld a = std::fabs(x);
    if (!(a > 0.5L) || a > 1) return std::acos(x);
    ld th, tl;
    asin_acos_t(a, th, tl);
    return std::atan2(th, x) + x * tl;
}
}  // namespace x87

// ---- overloads (preferred to the forwarders in quaternion.hpp for long double) ----
[[nodiscard]] inline long double sin(long double x) noexcept {
    return std::isfinite(x) ? x87::trig_eval(x, x87::trig::sin) : std::sin(x);
}
[[nodiscard]] inline long double cos(long double x) noexcept {
    return std::isfinite(x) ? x87::trig_eval(x, x87::trig::cos) : std::cos(x);
}
[[nodiscard]] inline long double tan(long double x) noexcept {
    return std::isfinite(x) ? x87::trig_eval(x, x87::trig::tan) : std::tan(x);
}
[[nodiscard]] inline long double sinh(long double x) noexcept { return x87::sinh(x); }
[[nodiscard]] inline long double expm1(long double x) noexcept { return x87::expm1(x); }
[[nodiscard]] inline long double pow(long double x, long double y) noexcept { return x87::pow(x, y); }
// pow(long double, double), pow(int, long double), ...: what std::pow computes for these (both
// arguments converted to long double); without this the forwarder would be the better match.
template <class A, class B>
    requires std::is_arithmetic_v<A> && std::is_arithmetic_v<B> && (!std::is_same_v<A, B>) &&
             (std::is_same_v<A, long double> || std::is_same_v<B, long double>)
[[nodiscard]] inline long double pow(A x, B y) noexcept {
    return x87::pow(static_cast<long double>(x), static_cast<long double>(y));
}
[[nodiscard]] inline long double asin(long double x) noexcept { return x87::asin(x); }
[[nodiscard]] inline long double acos(long double x) noexcept { return x87::acos(x); }
[[nodiscard]] inline long double tgamma(long double x) noexcept { return x87::tgamma(x); }
[[nodiscard]] inline long double lgamma(long double x) noexcept { return x87::lgamma(x); }

// ---- std::complex<long double> -----------------------------------------------------
// libc++ builds these on sinl, cosl and sinhl.  For finite arguments, the same formulas as
// libc++ with the functions above; non-finite arguments (whose special cases need no accurate
// sin or cos) go to std:: unchanged.
using cld = std::complex<long double>;
namespace x87 {
[[nodiscard]] inline bool finite(const std::complex<ld>& z) noexcept { return std::isfinite(z.real()) && std::isfinite(z.imag()); }
}  // namespace x87
[[nodiscard]] inline cld polar(long double rho, long double theta) noexcept {
    if (!std::isfinite(rho) || !std::isfinite(theta) || std::signbit(rho)) return std::polar(rho, theta);
    long double s, c;
    x87::sincos(theta, s, c);
    return {rho * c, rho * s};
}
[[nodiscard]] inline cld exp(const cld& z) noexcept {
    if (!x87::finite(z)) return std::exp(z);
    if (z.imag() == 0) return {std::exp(z.real()), z.imag()};
    const long double e = std::exp(z.real());
    long double s, c;
    x87::sincos(z.imag(), s, c);
    return {e * c, e * s};
}
[[nodiscard]] inline cld sinh(const cld& z) noexcept {
    if (!x87::finite(z)) return std::sinh(z);
    long double s, c;
    x87::sincos(z.imag(), s, c);
    return {lm::sinh(z.real()) * c, std::cosh(z.real()) * s};
}
[[nodiscard]] inline cld cosh(const cld& z) noexcept {
    if (!x87::finite(z)) return std::cosh(z);
    if (z.real() == 0 && z.imag() == 0) return {1.0L, z.imag()};
    long double s, c;
    x87::sincos(z.imag(), s, c);
    return {std::cosh(z.real()) * c, lm::sinh(z.real()) * s};
}
[[nodiscard]] inline cld sin(const cld& z) noexcept {  // -i sinh(i z)
    if (!x87::finite(z)) return std::sin(z);
    const cld w = lm::sinh(cld(-z.imag(), z.real()));
    return {w.imag(), -w.real()};
}
[[nodiscard]] inline cld cos(const cld& z) noexcept {  // cosh(i z)
    if (!x87::finite(z)) return std::cos(z);
    return lm::cosh(cld(-z.imag(), z.real()));
}
[[nodiscard]] inline cld sqrt(const cld& z) noexcept {
    if (!x87::finite(z)) return std::sqrt(z);
    return lm::polar(std::sqrt(std::abs(z)), std::arg(z) / 2.0L);
}
}  // namespace quat::detail::lm

#endif  // QUATERNION_DETAIL_X87_LIBM_HPP
