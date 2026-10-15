// The user turns contraction on before including the headers.
#include "common.hpp"
#if QFPC_MSVC
#pragma fp_contract(on)
#endif
QFPC_NOINLINE double qfpc::on_before(double a, double b, double c) { return a * b + c; }
#include <quaternion.hpp>
QFPC_NOINLINE double qfpc::on_after(double a, double b, double c) { return a * b + c; }
double qfpc_on_use(double x) { return quat::abs(quat::exp(quat::quaternion<double>(x, 1, 2, 3))); }
