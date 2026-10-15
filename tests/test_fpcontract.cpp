// The headers must not change the includer's floating-point contraction state.  MSVC with
// /fp:contract (or /fp:fast) contracts a * b + c into an fma; the headers compile their own
// code with #pragma fp_contract(off) there and restore the includer's state at their end
// (float_control(push)/(pop)).  Three translation units: this one (state from the command
// line), tests/fpcontract/off.cpp (#pragma fp_contract(off) before the include) and
// tests/fpcontract/on.cpp (#pragma fp_contract(on)).  In each, a * b + c defined after the
// include must give the bits of the same expression defined before it.  The inputs make the
// fused and unfused results differ: a * b = 1 - 2^-60 exactly, which rounds to 1, so a * b + c
// is 0 unfused and -2^-60 fused.
#include "fpcontract/common.hpp"
#include <cstdint>
#include <cstdio>
#include <cstring>
namespace qfpc {
QFPC_NOINLINE double cmd_before(double a, double b, double c) { return a * b + c; }
}
#include <quaternion.hpp>
namespace qfpc {
QFPC_NOINLINE double cmd_after(double a, double b, double c) { return a * b + c; }
}
double qfpc_off_use(double);
double qfpc_on_use(double);

static std::uint64_t bits(double x) { std::uint64_t u; std::memcpy(&u, &x, sizeof u); return u; }

int main() {
    volatile double va = 1 + 0x1p-30, vb = 1 - 0x1p-30, vc = -1;
    const double a = va, b = vb, c = vc;
    int fail = 0;
    struct { const char* state; double (*before)(double, double, double); double (*after)(double, double, double); } t[] = {
        {"command line", qfpc::cmd_before, qfpc::cmd_after},
        {"#pragma fp_contract(off)", qfpc::off_before, qfpc::off_after},
        {"#pragma fp_contract(on)", qfpc::on_before, qfpc::on_after},
    };
    bool fused[3];
    for (int i = 0; i < 3; ++i) {
        const double x = t[i].before(a, b, c), y = t[i].after(a, b, c);
        fused[i] = x != 0;
        std::printf("  %-26s before include: %s, after include: %s\n", t[i].state, x != 0 ? "fused" : "not fused",
                    y != 0 ? "fused" : "not fused");
        if (bits(x) != bits(y)) { std::printf("FAIL: the headers changed the contraction state (%s)\n", t[i].state); ++fail; }
    }
#if QFPC_MSVC
    // The test only means something if MSVC honours the states: with fp_contract(off) nothing
    // is fused; with /fp:contract and an fma instruction (x64 /arch:AVX2, ARM64) the
    // command-line state and fp_contract(on) are.
    if (fused[1]) { std::printf("FAIL: fp_contract(off) before the include still fused\n"); ++fail; }
#if defined(_M_FP_CONTRACT) && (defined(__AVX2__) || defined(_M_ARM64))
    if (!fused[0] || !fused[2]) { std::printf("FAIL: /fp:contract build did not fuse before the include\n"); ++fail; }
#endif
#endif
    (void)fused;
    std::printf("%s (%g, %g)\n", fail ? "FAILED" : "OK", qfpc_off_use(0.5), qfpc_on_use(0.5));
    return fail ? 1 : 0;
}
