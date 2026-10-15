// Shared by the fp-contract state test (tests/test_fpcontract.cpp): each translation unit
// defines f_before (above its #include <quaternion.hpp>) and f_after (below it), both
// returning a * b + c, under the same user pragma state.  The headers must leave that state
// as they found it, so the two must give the same bits.
#pragma once
#if defined(_MSC_VER)
#define QFPC_NOINLINE __declspec(noinline)
#else
#define QFPC_NOINLINE __attribute__((noinline))
#endif
// MSVC proper (not clang-cl): the compiler whose contraction the headers switch off.
#if defined(_MSC_VER) && !defined(__clang__)
#define QFPC_MSVC 1
#else
#define QFPC_MSVC 0
#endif
namespace qfpc {
double off_before(double, double, double);
double off_after(double, double, double);
double on_before(double, double, double);
double on_after(double, double, double);
}  // namespace qfpc
