// Unit tests for the GF(2^8) erasure code engine (MP-QKD-TP §V-I/§V-J)
//
// Compile (do not chain with run):
//   g++ -std=c++17 -O3 -mpclmul -Wall -Wextra -Iinclude src/*.cpp tests/test_main.cpp -o build/bin/qkdpp_tests.exe
// Then run:
//   ./build/bin/qkdpp_tests.exe

#include "qkdpp/gf2_8.hpp"
#include "qkdpp/vandermonde.hpp"
#include "qkdpp/erasure_code.hpp"

#include <cstdio>
#include <cstdint>
#include <vector>
#include <random>
#include <algorithm>
#include <stdexcept>

using namespace qkdpp;

static int g_tests_run    = 0;
static int g_tests_passed = 0;
static int g_tests_failed = 0;

static std::vector<void(*)()> g_test_funcs;
static std::vector<const char*> g_test_names;

#define TEST(name) \
    static void name(); \
    struct name##_Registerer { \
        name##_Registerer() { g_test_funcs.push_back(name); g_test_names.push_back(#name); } \
    }; \
    static name##_Registerer name##_reg; \
    static void name()

// Simple assertion helpers
static void check(bool cond, const char* msg, const char* file, int line) {
    if (!cond) {
        printf("  [FAIL] %s (%s:%d)\n", msg, file, line);
        g_tests_failed++;
    } else {
        g_tests_passed++;
    }
}

#define CHECK(cond, msg) check(cond, msg, __FILE__, __LINE__)
#define CHECK_EQ(a, b, msg) check((a) == (b), msg, __FILE__, __LINE__)

// ========================================================================
// GF(2^8) tests
// ========================================================================

TEST(test_gf_add_xor) {
    GF256 f;
    CHECK_EQ(f.add(0x53, 0xCA), (uint8_t)(0x53 ^ 0xCA), "add == XOR");
    CHECK_EQ(f.add(0x00, 0x00), (uint8_t)0, "0+0=0");
    CHECK_EQ(f.add(0xFF, 0xFF), (uint8_t)0, "a+a=0 (char 2)");
    CHECK_EQ(f.sub(0x53, 0xCA), (uint8_t)(0x53 ^ 0xCA), "sub == XOR");
}

TEST(test_gf_mul_table_vs_slow) {
    GF256 f;
    bool all_ok = true;
    // Test all 256x256 products
    for (int a = 0; a < 256; a++) {
        for (int b = 0; b < 256; b++) {
            uint8_t fast = f.mul(static_cast<uint8_t>(a), static_cast<uint8_t>(b));
            uint8_t slow = GF256::mul_slow(static_cast<uint8_t>(a), static_cast<uint8_t>(b));
            if (fast != slow) {
                all_ok = false;
                printf("    mismatch: mul(%d,%d) fast=%d slow=%d\n", a, b, fast, slow);
            }
        }
    }
    CHECK(all_ok, "mul table matches mul_slow for all 65536 inputs");
}

TEST(test_gf_mul_zero) {
    GF256 f;
    CHECK_EQ(f.mul(0, 5), (uint8_t)0, "0*5=0");
    CHECK_EQ(f.mul(5, 0), (uint8_t)0, "5*0=0");
    CHECK_EQ(f.mul(0, 0), (uint8_t)0, "0*0=0");
}

TEST(test_gf_mul_one) {
    GF256 f;
    CHECK_EQ(f.mul(1, 5), (uint8_t)5, "1*5=5");
    CHECK_EQ(f.mul(5, 1), (uint8_t)5, "5*1=5");
    CHECK_EQ(f.mul(1, 1), (uint8_t)1, "1*1=1");
}

TEST(test_gf_mul_commutativity) {
    GF256 f;
    bool ok = true;
    for (int a = 0; a < 256 && ok; a++) {
        for (int b = 0; b < 256; b++) {
            if (f.mul(a, b) != f.mul(b, a)) { ok = false; break; }
        }
    }
    CHECK(ok, "mul is commutative");
}

TEST(test_gf_mul_associativity) {
    GF256 f;
    bool ok = true;
    for (int a = 1; a < 256 && ok; a += 7) {
        for (int b = 1; b < 256 && ok; b += 11) {
            for (int c = 1; c < 256 && ok; c += 13) {
                uint8_t left  = f.mul(f.mul(a, b), c);
                uint8_t right = f.mul(a, f.mul(b, c));
                if (left != right) { ok = false; }
            }
        }
    }
    CHECK(ok, "mul is associative (sampled)");
}

TEST(test_gf_distributivity) {
    GF256 f;
    bool ok = true;
    for (int a = 1; a < 256 && ok; a += 5) {
        for (int b = 0; b < 256 && ok; b += 7) {
            for (int c = 0; c < 256 && ok; c += 11) {
                uint8_t left  = f.mul(a, f.add(b, c));
                uint8_t right = f.add(f.mul(a, b), f.mul(a, c));
                if (left != right) { ok = false; }
            }
        }
    }
    CHECK(ok, "mul distributes over add (sampled)");
}

TEST(test_gf_inverse_identity) {
    GF256 f;
    bool ok = true;
    for (int a = 1; a < 256; a++) {
        uint8_t inv = f.inverse(static_cast<uint8_t>(a));
        uint8_t prod = f.mul(static_cast<uint8_t>(a), inv);
        if (prod != 1) {
            ok = false;
            printf("    a=%d inv=%d a*inv=%d\n", a, inv, prod);
            break;
        }
    }
    CHECK(ok, "a * a^-1 = 1 for all nonzero a");
}

TEST(test_gf_inverse_double) {
    GF256 f;
    // inverse(inverse(a)) == a
    bool ok = true;
    for (int a = 1; a < 256; a++) {
        uint8_t inv1 = f.inverse(static_cast<uint8_t>(a));
        uint8_t inv2 = f.inverse(inv1);
        if (inv2 != a) { ok = false; break; }
    }
    CHECK(ok, "inverse(inverse(a)) == a");
}

TEST(test_gf_div) {
    GF256 f;
    // a / b == a * b^-1
    bool ok = true;
    for (int a = 1; a < 256 && ok; a += 3) {
        for (int b = 1; b < 256 && ok; b += 5) {
            uint8_t quot = f.div(static_cast<uint8_t>(a), static_cast<uint8_t>(b));
            uint8_t expected = f.mul(static_cast<uint8_t>(a),
                                      f.inverse(static_cast<uint8_t>(b)));
            if (quot != expected) { ok = false; }
        }
    }
    CHECK(ok, "div(a,b) == a * b^-1 (sampled)");
}

TEST(test_gf_pow) {
    GF256 f;
    // 3^255 = 1 (Fermat-like, since the multiplicative group has order 255)
    CHECK_EQ(f.pow(3, 255), (uint8_t)1, "alpha^255 = 1");
    CHECK_EQ(f.pow(3, 0), (uint8_t)1, "alpha^0 = 1");
    CHECK_EQ(f.pow(3, 1), (uint8_t)3, "alpha^1 = 3");
    // alpha^2 = 3*3 via mul_slow = 5
    CHECK_EQ(f.pow(3, 2), GF256::mul_slow(3, 3), "alpha^2 correct");
    // a^256 = a for any nonzero a (a^255 = 1, so a^256 = a)
    bool ok = true;
    for (int a = 1; a < 256 && ok; a++) {
        if (f.pow(static_cast<uint8_t>(a), 256) != a) ok = false;
    }
    CHECK(ok, "a^256 = a for all nonzero a");
}

TEST(test_gf_log_antilog) {
    GF256 f;
    // alpha^log(x) = x for x in [1,255]
    bool ok = true;
    for (int x = 1; x < 256; x++) {
        uint8_t lg = f.log_val(static_cast<uint8_t>(x));
        uint8_t recovered = f.antilog(lg);
        if (recovered != x) { ok = false; break; }
    }
    CHECK(ok, "alpha^log(x) = x for all nonzero x");
}

// ========================================================================
// Vandermonde matrix tests
// ========================================================================

TEST(test_vandermonde_shape) {
    GF256 f;
    auto V = Vandermonde::build(f, 4, 2);
    CHECK_EQ(static_cast<int>(V.size()), 4, "4 rows");
    CHECK_EQ(static_cast<int>(V[0].size()), 2, "2 cols");
}

TEST(test_vandermonde_first_row_all_ones) {
    GF256 f;
    auto V = Vandermonde::build(f, 4, 4);
    // Row 0: [1, 1, 1, 1] since alpha^(0*j) = 1
    bool ok = true;
    for (int j = 0; j < 4; j++) {
        if (V[0][j] != 1) ok = false;
    }
    CHECK(ok, "first row is all 1s");
}

TEST(test_vandermonde_first_col_all_ones) {
    GF256 f;
    auto V = Vandermonde::build(f, 4, 3);
    // Col 0: alpha^(i*0) = 1 for all i
    bool ok = true;
    for (int i = 0; i < 4; i++) {
        if (V[i][0] != 1) ok = false;
    }
    CHECK(ok, "first column is all 1s");
}

TEST(test_vandermonde_independence_K4) {
    GF256 f;
    // K' x K' submatrix (first K' rows) must be linearly independent
    auto V = Vandermonde::build(f, 4, 4);
    CHECK(Vandermonde::rows_linearly_independent(f, V), "4x4 Vandermonde rows independent");
}

TEST(test_vandermonde_any_submatrix_independent) {
    GF256 f;
    // Any K' x K' submatrix of a Vandermonde must be non-singular (MDS)
    // Test with K=6, K'=4 and check several row subsets.
    auto V = Vandermonde::build(f, 6, 4);

    // Check {0,1,2,3}, {0,2,4,5}, {1,3,4,5}, {0,1,4,5}
    std::vector<std::vector<int>> subsets = {
        {0,1,2,3}, {0,2,4,5}, {1,3,4,5}, {0,1,4,5}
    };
    bool all_ok = true;
    for (const auto& s : subsets) {
        std::vector<std::vector<uint8_t>> sub(4, std::vector<uint8_t>(4, 0));
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                sub[i][j] = V[s[i]][j];
            }
        }
        if (!Vandermonde::rows_linearly_independent(f, sub)) {
            all_ok = false;
            printf("    subset {");
            for (int x : s) printf("%d ", x);
            printf("} is singular\n");
        }
    }
    CHECK(all_ok, "any K'xK' submatrix of Vandermonde is non-singular");
}

// ========================================================================
// Erasure code tests
// ========================================================================

// Helper: verify decode round-trip
static bool verify_round_trip(
    int K, int K_prime,
    const std::vector<std::vector<uint8_t>>& original,
    const std::vector<int>& surviving_indices) {

    auto encoded = ErasureCode::encode(K, K_prime, original);

    // Extract surviving rows
    std::vector<std::vector<uint8_t>> surviving;
    for (int idx : surviving_indices) {
        surviving.push_back(encoded[static_cast<size_t>(idx)]);
    }

    auto recovered = ErasureCode::decode(K, K_prime, surviving_indices, surviving);

    if (recovered.size() != original.size()) return false;
    for (size_t i = 0; i < original.size(); i++) {
        if (recovered[i] != original[i]) return false;
    }
    return true;
}

TEST(test_erasure_encode_manual_K4_K2) {
    // K=4, K'=2, verify encoded rows against hand computation
    GF256 f;
    // data: 2 rows of length 1
    // d0 = [0x01], d1 = [0x02]
    // Vandermonde V is 4x2:
    //   row0: [1, 1]
    //   row1: [1, alpha] = [1, 3]
    //   row2: [1, alpha^2] = [1, 5]   (3*3=5 via mul_slow)
    //   row3: [1, alpha^3] = [1, 15]  (5*3=15 via mul_slow)
    // encoded[i] = V[i][0]*d0 + V[i][1]*d1 = 1*0x01 + alpha^i*0x02
    //   e0 = 1 ^ (1*2) = 1 ^ 2 = 3
    //   e1 = 1 ^ (3*2) = 1 ^ 6 = 7
    //   e2 = 1 ^ (5*2) = 1 ^ 10 = 11
    //   e3 = 1 ^ (15*2) = 1 ^ 30 = 31

    std::vector<std::vector<uint8_t>> data = {{0x01}, {0x02}};
    auto encoded = ErasureCode::encode(4, 2, data);

    CHECK_EQ(static_cast<int>(encoded.size()), 4, "4 encoded rows");
    CHECK_EQ(static_cast<int>(encoded[0].size()), 1, "row length 1");

    // Verify expected values
    CHECK_EQ(encoded[0][0], (uint8_t)(0x01 ^ 0x02), "e0 = 3");
    uint8_t e1 = 0x01 ^ GF256::mul_slow(3, 0x02);
    CHECK_EQ(encoded[1][0], e1, "e1 = 7");
    uint8_t e2 = 0x01 ^ GF256::mul_slow(5, 0x02);
    CHECK_EQ(encoded[2][0], e2, "e2 = 11");
    uint8_t e3 = 0x01 ^ GF256::mul_slow(15, 0x02);
    CHECK_EQ(encoded[3][0], e3, "e3 = 31");
}

TEST(test_erasure_decode_no_loss_K4_K2) {
    // No loss: use first K' rows
    std::vector<std::vector<uint8_t>> data = {
        {0x42, 0x01, 0xFF},
        {0x80, 0x7F, 0x00}
    };
    bool ok = verify_round_trip(4, 2, data, {0, 1});
    CHECK(ok, "no loss: decode(0,1) == original");
}

TEST(test_erasure_decode_one_loss_K4_K2) {
    // Lose row 0, use rows {1, 2}
    std::vector<std::vector<uint8_t>> data = {
        {0x42, 0x01, 0xFF, 0xAA},
        {0x80, 0x7F, 0x00, 0x55}
    };
    bool ok = verify_round_trip(4, 2, data, {1, 2});
    CHECK(ok, "1 loss: decode(1,2) == original");
}

TEST(test_erasure_decode_max_loss_K4_K2) {
    // K - K' = 2 losses (maximum), use rows {2, 3}
    std::vector<std::vector<uint8_t>> data = {
        {0x42, 0x01, 0xFF, 0xAA, 0xBB},
        {0x80, 0x7F, 0x00, 0x55, 0xCC}
    };
    bool ok = verify_round_trip(4, 2, data, {2, 3});
    CHECK(ok, "max loss: decode(2,3) == original");
}

TEST(test_erasure_decode_all_subsets_K4_K2) {
    // Test all C(4,2) = 6 combinations of surviving rows
    std::vector<std::vector<uint8_t>> data = {
        {0x12, 0x34, 0x56, 0x78, 0x9A, 0xBC, 0xDE, 0xF0},
        {0xFE, 0xDC, 0xBA, 0x98, 0x76, 0x54, 0x32, 0x10}
    };
    int subsets[][2] = {
        {0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}
    };
    bool all_ok = true;
    for (auto& s : subsets) {
        if (!verify_round_trip(4, 2, data, {s[0], s[1]})) {
            all_ok = false;
            printf("    failed for subset {%d,%d}\n", s[0], s[1]);
        }
    }
    CHECK(all_ok, "all 6 subsets of 2-out-of-4 decode correctly");
}

TEST(test_erasure_decode_random_data_K4_K2) {
    // Random data, multiple trials
    std::mt19937 rng(12345);
    std::vector<std::vector<uint8_t>> data(2, std::vector<uint8_t>(16, 0));
    bool all_ok = true;
    for (int trial = 0; trial < 50; trial++) {
        for (auto& row : data) {
            for (auto& b : row) {
                b = static_cast<uint8_t>(rng());
            }
        }
        // Randomly pick 2 of 4 rows
        int s0 = rng() % 4;
        int s1;
        do { s1 = rng() % 4; } while (s1 == s0);
        if (s0 > s1) std::swap(s0, s1);

        if (!verify_round_trip(4, 2, data, {s0, s1})) {
            all_ok = false;
            printf("    trial %d failed for subset {%d,%d}\n", trial, s0, s1);
        }
    }
    CHECK(all_ok, "50 random trials with K=4,K'=2 all decode correctly");
}

TEST(test_erasure_larger_K8_K4) {
    // K=8, K'=4, matches paper's "half rows lost" scenario
    std::mt19937 rng(42);
    std::vector<std::vector<uint8_t>> data(4, std::vector<uint8_t>(32, 0));
    for (auto& row : data) {
        for (auto& b : row) {
            b = static_cast<uint8_t>(rng());
        }
    }

    // Test losing 0, 1, 2, 3, 4 (= K - K') rows
    bool all_ok = true;
    for (int lost = 0; lost <= 4 && all_ok; lost++) {
        // Pick surviving indices: take rows [lost .. lost+3] (mod 8)
        std::vector<int> surviving;
        for (int i = 0; i < 4; i++) {
            surviving.push_back((lost + i) % 8);
        }
        // Ensure distinct
        std::sort(surviving.begin(), surviving.end());
        auto last = std::unique(surviving.begin(), surviving.end());
        if (static_cast<int>(std::distance(surviving.begin(), last)) != 4) {
            continue;
        }
        if (!verify_round_trip(8, 4, data, surviving)) {
            all_ok = false;
            printf("    failed for lost=%d\n", lost);
        }
    }
    CHECK(all_ok, "K=8,K'=4: 0..4 losses all decode correctly");
}

TEST(test_erasure_decode_max_loss_K8_K4) {
    // K=8, K'=4: lose exactly 4 (= K-K') rows, use last 4 rows {4,5,6,7}
    std::mt19937 rng(999);
    std::vector<std::vector<uint8_t>> data(4, std::vector<uint8_t>(64, 0));
    for (auto& row : data) {
        for (auto& b : row) {
            b = static_cast<uint8_t>(rng());
        }
    }
    bool ok = verify_round_trip(8, 4, data, {4, 5, 6, 7});
    CHECK(ok, "K=8,K'=4: max loss (rows 4-7) decodes correctly");
}

TEST(test_erasure_multi_byte_row) {
    // Rows of length 256, all values
    std::vector<std::vector<uint8_t>> data(2, std::vector<uint8_t>(256));
    for (int b = 0; b < 256; b++) {
        data[0][b] = static_cast<uint8_t>(b);
        data[1][b] = static_cast<uint8_t>(255 - b);
    }
    bool ok = verify_round_trip(4, 2, data, {2, 3});
    CHECK(ok, "256-byte rows, max loss, decodes correctly");
}

TEST(test_erasure_single_byte_rows) {
    // Edge: rows of length 1
    std::vector<std::vector<uint8_t>> data = {{0x5A}, {0xA5}};
    bool ok = verify_round_trip(4, 2, data, {1, 3});
    CHECK(ok, "1-byte rows, scattered loss, decodes correctly");
}

TEST(test_erasure_empty_rows) {
    // Edge: rows of length 0
    std::vector<std::vector<uint8_t>> data = {{}, {}};
    bool ok = verify_round_trip(4, 2, data, {1, 2});
    CHECK(ok, "0-byte rows decode to 0-byte rows");
}

// --- Dimension limits (GF(2^8)* has order 255) ---------------------------

TEST(test_vandermonde_cols_upper_bound) {
    GF256 f;
    // 255 columns is the maximum distinct-power set
    bool ok255 = true;
    try { (void)Vandermonde::build(f, 4, 255); } catch (...) { ok255 = false; }
    CHECK(ok255, "cols == 255 is accepted");

    // 256 columns would duplicate column 0 -> must be rejected
    bool threw = false;
    try { (void)Vandermonde::build(f, 4, 256); }
    catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw, "cols == 256 is rejected (period 255)");
}

TEST(test_vandermonde_rows_upper_bound) {
    GF256 f;
    bool ok255 = true;
    try { (void)Vandermonde::build(f, 255, 4); } catch (...) { ok255 = false; }
    CHECK(ok255, "rows == 255 is accepted");

    // row 255 would repeat row 0
    bool threw = false;
    try { (void)Vandermonde::build(f, 256, 4); }
    catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw, "rows == 256 is rejected (period 255)");

    // large dimensions must not overflow the exponent computation
    bool ok_big = true;
    try {
        auto M = Vandermonde::build(f, 255, 255);
        ok_big = (M.size() == 255 && M[0].size() == 255);
        // row 1 = alpha^0 .. alpha^254 -> all distinct non-zero values
        if (ok_big) {
            std::vector<bool> seen(256, false);
            for (int j = 0; j < 255; ++j) {
                uint8_t v = M[1][static_cast<size_t>(j)];
                if (v == 0 || seen[v]) ok_big = false;
                seen[v] = true;
            }
        }
    } catch (...) { ok_big = false; }
    CHECK(ok_big, "255x255 Vandermonde builds with 255 distinct values in row 1");
}

TEST(test_vandermonde_ragged_matrix_rejected) {
    GF256 f;
    std::vector<std::vector<uint8_t>> ragged = {{1, 2, 3}, {4, 5}};
    bool threw = false;
    try { (void)Vandermonde::rows_linearly_independent(f, ragged); }
    catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw, "ragged matrix is rejected instead of reading out of bounds");
}

TEST(test_erasure_dimension_limits) {
    std::vector<std::vector<uint8_t>> data(2, std::vector<uint8_t>(1, 0x11));

    bool threw_k = false;
    try { (void)ErasureCode::encode(256, 2, data); }
    catch (const std::invalid_argument&) { threw_k = true; }
    CHECK(threw_k, "encode rejects K > 255");

    bool threw_kp = false;
    try { (void)ErasureCode::encode(300, 256, std::vector<std::vector<uint8_t>>(256, std::vector<uint8_t>(1, 0))); }
    catch (const std::invalid_argument&) { threw_kp = true; }
    CHECK(threw_kp, "encode rejects K' > 255");

    bool threw_dec = false;
    try { (void)ErasureCode::decode(256, 2, {0, 1}, data); }
    catch (const std::invalid_argument&) { threw_dec = true; }
    CHECK(threw_dec, "decode rejects K > 255");

    // K == 255 is the largest admissible configuration
    bool ok_max = true;
    try {
        std::vector<std::vector<uint8_t>> d(2, std::vector<uint8_t>(3, 0x7E));
        auto enc = ErasureCode::encode(255, 2, d);
        auto dec = ErasureCode::decode(255, 2, {253, 254}, {enc[253], enc[254]});
        ok_max = (dec == d);
    } catch (...) { ok_max = false; }
    CHECK(ok_max, "K == 255 round-trips through encode/decode");
}

// ========================================================================
// Main
// ========================================================================

int main() {
    printf("=== libqkd-erasure unit tests ===\n\n");

    for (size_t i = 0; i < g_test_funcs.size(); i++) {
        printf("--- %s ---\n", g_test_names[i]);
        g_tests_run++;
        g_test_funcs[i]();
    }

    printf("\n=== Summary ===\n");
    printf("Tests run:    %d\n", g_tests_run);
    printf("Checks passed: %d\n", g_tests_passed);
    printf("Checks failed: %d\n", g_tests_failed);

    if (g_tests_failed == 0) {
        printf("\nALL TESTS PASSED\n");
        return 0;
    } else {
        printf("\nSOME TESTS FAILED\n");
        return 1;
    }
}
