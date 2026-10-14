// Prints the bit-identity hashes ("raw" bits, or "canonical" NaNs), one "type op hash" per line.
#include "ops.hpp"
#include <cstdio>
#include <cstring>
template <class T>
void dump(const char* t, bool c) {
    for (const auto& [name, h] : bitident::run<T>(c)) std::printf("%s %s %016llx\n", t, name.c_str(), static_cast<unsigned long long>(h));
}
int main(int argc, char** argv) {
    const bool c = argc > 1 && std::strcmp(argv[1], "canonical") == 0;
    dump<float>("float", c);
    dump<double>("double", c);
    dump<long double>("long_double", c);
}
