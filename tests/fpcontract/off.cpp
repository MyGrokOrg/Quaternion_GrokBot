// The user turns contraction off before including the headers.
#include "common.hpp"
#if QFPC_MSVC
#pragma fp_contract(off)
#endif
QFPC_NOINLINE double qfpc::off_before(double a, double b, double c) { return a * b + c; }
#include <quaternion.hpp>
QFPC_NOINLINE double qfpc::off_after(double a, double b, double c) { return a * b + c; }
// Keep the headers' code instantiated in this TU too.
double qfpc_off_use(double x) { return quat::abs(quat::exp(quat::quaternion<double>(x, 1, 2, 3))); }
