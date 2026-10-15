// Runtime results of the transcendental-free operations must stay bit-identical to the
// header before the constexpr work (golden hashes in bitident/golden.inc; NaNs canonical).
// With FMA (e.g. -march=haswell) the contraction pattern depends on the compiler and the
// optimization level, so no golden hash is stable there; those builds are checked instead by
// the CI job "bit-identity (haswell)", which diffs the raw bits of the old and new header
// compiled with the same compiler and flags (tests/bitident/dump.cpp).
#include "test_framework.hpp"
#include "bitident/ops.hpp"
#include <cstring>
#include <limits>
#include <string>

struct golden { const char* type; const char* op; std::uint64_t hash; };
static const golden G[] = {
#include "bitident/golden.inc"
};

// abs and normalize call the C library's hypot.  The goldens are glibc's; the Universal CRT
// (MSVC, clang-cl) rounds hypot differently, so Windows has its own hashes for those two
// ("<type>@ucrt" entries in golden.inc).
[[maybe_unused]] static bool uses_libm_hypot(const std::string& op) { return op == "abs" || op == "normalize"; }

template <class T>
void check(const char* tname) {
    // Where long double is binary64 (MSVC), its results must equal double's bit for bit.
    if (std::numeric_limits<T>::digits == 53 && std::strcmp(tname, "long_double") == 0) tname = "double";
    for (const auto& [name, h] : bitident::run<T>(true)) {
        bool found = false;
        std::string type = tname;
#if defined(_WIN32)
        if (uses_libm_hypot(name)) type += "@ucrt";
#endif
        for (const auto& g : G) {
            if (type != g.type || name != g.op) continue;
            found = true;
            const std::uint64_t want = g.hash;
            if (h != want) std::printf("  %s %s: %016llx, golden %016llx\n", tname, name.c_str(),
                                       static_cast<unsigned long long>(h), static_cast<unsigned long long>(want));
            CHECK(h == want);
        }
        CHECK(found);
    }
}
TEST(runtime_bit_identity) {
#if defined(__FMA__)
    std::printf("  FMA build: golden hashes skipped (see the bit-identity CI job)\n");
    return;
#endif
    check<float>("float");
    check<double>("double");
    check<long double>("long_double");
}
TEST_MAIN
