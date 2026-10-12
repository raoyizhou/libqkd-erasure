#ifndef QKDPP_GF2_8_HPP
#define QKDPP_GF2_8_HPP

#include <cstdint>
#include <array>

namespace qkdpp {

// GF(2^8) field arithmetic using the AES polynomial 0x11B
// (x^8 + x^4 + x^3 + x + 1).
// Uses log/antilog (Zech) tables with generator alpha = 0x03 for O(1) multiply.
class GF256 {
public:
    static constexpr uint16_t POLY = 0x11B;   // x^8 + x^4 + x^3 + x + 1
    static constexpr uint8_t  ALPHA = 0x03;    // primitive element
    // Order of the multiplicative group GF(2^8)* : the number of distinct
    // non-zero elements, and therefore the period of alpha^i.  Powers that
    // reach this value wrap back to alpha^0 = 1.
    static constexpr int      ORDER = 255;

    GF256();   // builds log/antilog tables

    // Field operations
    uint8_t add(uint8_t a, uint8_t b) const;       // XOR
    uint8_t sub(uint8_t a, uint8_t b) const;       // XOR (same as add)
    uint8_t mul(uint8_t a, uint8_t b) const;       // via log/antilog
    uint8_t div(uint8_t a, uint8_t b) const;       // a / b
    uint8_t pow(uint8_t a, int e) const;           // a^e
    uint8_t inverse(uint8_t a) const;              // multiplicative inverse

    // Direct table access (for testing / Vandermonde construction)
    uint8_t antilog(int i) const;   // alpha^i, i in [0,254] wraps to [0,255)
    uint8_t log_val(uint8_t x) const;  // log_alpha(x), x != 0

    // Russian-peasant multiply without table (for verification / bootstrap)
    static uint8_t mul_slow(uint8_t a, uint8_t b);

    // Build antilog table entry: alpha^i mod polynomial
    static uint8_t gen_pow(int i);

private:
    // log[x] = discrete log base alpha of x; 0 if x == 0 (sentinel)
    // antilog_table[i] = alpha^i for i in [0,255]; wraps around
    std::array<uint8_t, 256> log_table{};
    std::array<uint8_t, 256> antilog_table{};

    void init_tables();
};

} // namespace qkdpp

#endif // QKDPP_GF2_8_HPP
