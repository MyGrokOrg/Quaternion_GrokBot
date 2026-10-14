// Runtime results of the transcendental-free operations must stay bit-identical to the
// header before the constexpr work (golden hashes in bitident/golden.inc; NaNs canonical).
// With FMA (e.g. -march=haswell) the contraction pattern depends on the compiler and the
// optimization level, so no golden hash is stable there; those builds are checked instead by
// the CI job "bit-identity (haswell)", which diffs the raw bits of the old and new header
// compiled with the same compiler and flags (tests/bitident/dump.cpp).
#include "test_framework.hpp"
#include "bitident/ops.hpp"
#include <cstring>

struct golden { const char* type; const char* op; std::uint64_t hash; };
static const golden G[] = {
#include "bitident/golden.inc"
};

template <class T>
void check(const char* tname) {
    for (const auto& [name, h] : bitident::run<T>(true)) {
        bool found = false;
        for (const auto& g : G) {
            if (std::strcmp(g.type, tname) != 0 || name != g.op) continue;
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
