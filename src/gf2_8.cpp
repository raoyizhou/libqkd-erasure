#include "qkdpp/gf2_8.hpp"

namespace qkdpp {

// --- Static helpers -------------------------------------------------------

// Russian-peasant multiply in GF(2^8) with reduction by 0x11B.
// Used both as the bootstrap for table construction and as a slow
// independent implementation for testing the fast table path.
uint8_t GF256::mul_slow(uint8_t a, uint8_t b) {
    uint8_t p = 0;
    uint8_t acopy = a;
    uint8_t bcopy = b;
    for (int i = 0; i < 8; i++) {
        if (bcopy & 1) {
            p ^= acopy;
        }
        bool hi = (acopy & 0x80) != 0;
        acopy <<= 1;
        if (hi) {
            acopy ^= 0x1B; // low 8 bits of polynomial 0x11B
        }
        bcopy >>= 1;
    }
    return p;
}

// Compute alpha^i mod 0x11B (alpha = 0x03)
uint8_t GF256::gen_pow(int i) {
    i = ((i % 255) + 255) % 255; // wrap to [0,254]
    uint8_t result = 1;
    uint8_t base = ALPHA; // 0x03
    int e = i;
    while (e > 0) {
        if (e & 1) {
            result = mul_slow(result, base);
        }
        base = mul_slow(base, base);
        e >>= 1;
    }
    return result;
}

// --- Table initialization --------------------------------------------------

void GF256::init_tables() {
    // antilog[i] = alpha^i for i = 0..254; antilog[255] and the rest wrap.
    // log[x] = i such that alpha^i = x, for x in [1,255]; log[0] = 0 (sentinel)
    uint8_t x = 1; // alpha^0 = 1
    for (int i = 0; i < 255; i++) {
        antilog_table[i] = x;
        log_table[x] = static_cast<uint8_t>(i);
        x = mul_slow(x, ALPHA); // x = x * alpha
    }
    // For convenience, fill the wraparound so antilog[i] works for any i % 255
    for (int i = 255; i < 256; i++) {
        antilog_table[i] = antilog_table[i - 255]; // alpha^255 = alpha^0 = 1
    }
}

GF256::GF256() {
    init_tables();
}

// --- Field operations -----------------------------------------------------

uint8_t GF256::add(uint8_t a, uint8_t b) const {
    return a ^ b;
}

uint8_t GF256::sub(uint8_t a, uint8_t b) const {
    return a ^ b; // addition == subtraction in GF(2^8)
}

uint8_t GF256::mul(uint8_t a, uint8_t b) const {
    if (a == 0 || b == 0) return 0;
    int la = log_table[a];
    int lb = log_table[b];
    int sum = la + lb;
    // alpha^(la+lb) mod (2^8 - 1)
    if (sum >= 255) sum -= 255;
    return antilog_table[sum];
}

uint8_t GF256::div(uint8_t a, uint8_t b) const {
    if (a == 0) return 0;
    if (b == 0) return 0; // undefined; return 0 defensively
    int la = log_table[a];
    int lb = log_table[b];
    int diff = la - lb;
    if (diff < 0) diff += 255;
    return antilog_table[diff];
}

uint8_t GF256::pow(uint8_t a, int e) const {
    if (a == 0) {
        return (e == 0) ? 1 : 0;
    }
    // alpha^(log(a) * e) mod 255
    int la = log_table[a];
    // (la * e) mod 255
    long long prod = static_cast<long long>(la) * e;
    int r = static_cast<int>(prod % 255);
    if (r < 0) r += 255;
    return antilog_table[r];
}

uint8_t GF256::inverse(uint8_t a) const {
    if (a == 0) return 0; // 0 has no inverse; return 0 defensively
    int la = log_table[a];
    // a^-1 = alpha^(255 - la)
    int inv = 255 - la;
    if (inv == 255) inv = 0; // a == 1 -> inverse is 1
    return antilog_table[inv];
}

uint8_t GF256::antilog(int i) const {
    i = ((i % 255) + 255) % 255;
    return antilog_table[i];
}

uint8_t GF256::log_val(uint8_t x) const {
    return log_table[x];
}

} // namespace qkdpp
