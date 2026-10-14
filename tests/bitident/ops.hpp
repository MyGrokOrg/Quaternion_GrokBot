// Runtime bit-identity corpus for the transcendental-free operations (arithmetic incl.
// Annex G division, norm, conj, inverse, dot, cross, abs, normalize, proj, integer pow,
// comparisons).  Every result is hashed (FNV-1a over the value bits; in "canonical" mode
// every NaN hashes as one token, because NaN sign/payload may legitimately differ between
// compilers that order commutative operands differently).  Used by test_bitident (golden
// hashes from the pre-constexpr header) and by bitident_dump (old-vs-new header diff in CI).
// Only uses API that existed before the constexpr work, so it compiles against both headers.
#ifndef QUATERNION_BITIDENT_OPS_HPP
#define QUATERNION_BITIDENT_OPS_HPP
#include <quaternion.hpp>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

namespace bitident {

struct hasher {
    bool canonical;
    std::uint64_t h = 1469598103934665603ULL;
    void byte(unsigned char b) { h = (h ^ b) * 1099511628211ULL; }
    template <class T>
    void val(T x) {
        if (canonical && std::isnan(x)) { byte(0xA5); return; }
        unsigned char b[sizeof(T)] = {};
        std::memcpy(b, &x, sizeof(T));
        // x87 long double: only the 10 value bytes are significant.
        const std::size_t n = (std::numeric_limits<T>::digits == 64) ? 10 : sizeof(T);
        for (std::size_t i = 0; i < n; ++i) byte(b[i]);
    }
    template <class T>
    void q(const quat::quaternion<T>& v) { val(v.w()); val(v.x()); val(v.y()); val(v.z()); }
};

// Deterministic inputs: a 64-bit LCG drives sign, mantissa and a wide exponent range
// (including subnormal and near-overflow components), plus special values.
template <class T>
std::vector<quat::quaternion<T>> inputs(std::uint64_t seed, int n) {
    using L = std::numeric_limits<T>;
    std::uint64_t s = seed;
    const auto next = [&s] { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return s >> 11; };
    const auto comp = [&]() -> T {
        const std::uint64_t r = next();
        const int cls = static_cast<int>(r % 16);
        const T m = static_cast<T>(static_cast<double>(next() % (1ULL << 52)) / 4503599627370496.0 + 1.0);
        const T sg = (r >> 8) & 1 ? T(-1) : T(1);
        int e;
        if (cls < 8) e = static_cast<int>(next() % 40) - 20;                  // moderate
        else if (cls < 11) e = static_cast<int>(next() % static_cast<std::uint64_t>(2 * L::max_exponent)) - L::max_exponent;  // wide
        else if (cls < 12) e = L::max_exponent - 1 - static_cast<int>(next() % 4);  // near overflow
        else if (cls < 13) e = L::min_exponent - static_cast<int>(next() % static_cast<std::uint64_t>(L::digits + 2));  // subnormal
        else if (cls < 15) return sg * T(0);
        else return sg * static_cast<T>(static_cast<int>(next() % 7));  // small integers
        return sg * std::ldexp(m, e);
    };
    std::vector<quat::quaternion<T>> v;
    for (int i = 0; i < n; ++i) {
        const T a = comp(), b = comp(), c = comp(), d = comp();  // sequenced (argument order is unspecified)
        v.emplace_back(a, b, c, d);
    }
    const T inf = L::infinity(), nan = L::quiet_NaN();
    const T sp[] = {T(0), -T(0), T(1), T(-2.5), inf, -inf, nan, L::denorm_min(), L::max(), L::min()};
    for (T a : sp)
        for (T b : sp) v.emplace_back(a, b, T(0.5), -b);
    return v;
}

// One hash per operation name.
template <class T>
std::vector<std::pair<std::string, std::uint64_t>> run(bool canonical) {
    using Q = quat::quaternion<T>;
    volatile std::uint64_t seed = 0x5eed1234ULL;  // keep the inputs out of the constant folder
    const auto P = inputs<T>(seed, 1500), R = inputs<T>(seed ^ 0x9e3779b97f4a7c15ULL, 1500);
    std::vector<std::pair<std::string, std::uint64_t>> out;
    const auto op = [&](const char* name, auto f) {
        hasher h{canonical};
        for (std::size_t i = 0; i < P.size(); ++i) f(h, P[i], R[i]);
        out.emplace_back(name, h.h);
    };
    op("add", [](hasher& h, const Q& p, const Q& r) { h.q(p + r); h.q(p - r); h.q(p + r.w()); h.q(r.w() - p); });
    op("mul", [](hasher& h, const Q& p, const Q& r) { h.q(p * r); h.q(r * p); });
    op("mul_scalar", [](hasher& h, const Q& p, const Q& r) { h.q(p * r.x()); h.q(r.y() * p); h.q(p / r.z()); });
    op("div", [](hasher& h, const Q& p, const Q& r) { h.q(p / r); h.q(r / p); });
    op("scalar_div_q", [](hasher& h, const Q& p, const Q& r) { h.q(r.w() / p); });
    op("complex", [](hasher& h, const Q& p, const Q& r) {
        const std::complex<T> z(r.w(), r.x());
        h.q(p * z); h.q(z * p); h.q(p / z); h.q(z / p); h.q(p + z);
    });
    op("compound", [](hasher& h, const Q& p, const Q& r) { Q t = p; t *= r; t /= p; t += r; t -= p; h.q(t); });
    op("norm_dot", [](hasher& h, const Q& p, const Q& r) { h.val(norm(p)); h.val(dot(p, r)); });
    op("conj_cross", [](hasher& h, const Q& p, const Q& r) { h.q(conj(p)); h.q(cross(p, r)); h.q(-p); h.q(+p); });
    op("inverse", [](hasher& h, const Q& p, const Q&) { h.q(inverse(p)); });
    op("abs", [](hasher& h, const Q& p, const Q&) { h.val(abs(p)); });
    op("normalize", [](hasher& h, const Q& p, const Q&) { h.q(normalize(p)); });
    op("proj", [](hasher& h, const Q& p, const Q&) { h.q(proj(p)); });
    op("pow_int", [](hasher& h, const Q& p, const Q&) {
        for (int n = -3; n <= 5; ++n) h.q(pow(p, n));
    });
    op("compare", [](hasher& h, const Q& p, const Q& r) {
        h.byte(static_cast<unsigned char>((p == r) | ((p != r) << 1) | ((p == p.w()) << 2) | ((p == std::complex<T>(p.w(), p.x())) << 3)));
    });
    return out;
}

}  // namespace bitident
#endif
