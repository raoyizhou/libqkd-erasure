#include "qkdpp/erasure_code.hpp"
#include "qkdpp/gf2_8.hpp"
#include "qkdpp/vandermonde.hpp"
#include <stdexcept>
#include <algorithm>

namespace qkdpp {

// Upper bound for both K and K': the multiplicative group of GF(2^8) has
// order 255, so row i of the Vandermonde matrix repeats row (i - 255).
// Beyond that the matrix loses its MDS property and arbitrary K'xK'
// submatrices may be singular.
static int max_dimension() { return GF256::ORDER; }   // 255

std::vector<std::vector<uint8_t>> ErasureCode::encode(
    int K, int K_prime,
    const std::vector<std::vector<uint8_t>>& data) {

    if (K <= 0 || K_prime <= 0 || K < K_prime) {
        throw std::invalid_argument("ErasureCode::encode: need 1 <= K' <= K");
    }
    if (K_prime > max_dimension() || K > max_dimension()) {
        throw std::invalid_argument(
            "ErasureCode::encode: K and K' must be <= 255 for GF(2^8)");
    }
    if (static_cast<int>(data.size()) != K_prime) {
        throw std::invalid_argument(
            "ErasureCode::encode: data must have K' rows");
    }

    // Determine the row length (must be same for all data rows)
    size_t row_len = 0;
    if (!data.empty()) {
        row_len = data[0].size();
        for (const auto& row : data) {
            if (row.size() != row_len) {
                throw std::invalid_argument(
                    "ErasureCode::encode: all data rows must have same length");
            }
        }
    }

    GF256 field;
    auto V = Vandermonde::build(field, K, K_prime);

    // encoded = V * data, where V is K x K' and data is K' x row_len
    // Result: K x row_len
    std::vector<std::vector<uint8_t>> encoded(
        static_cast<size_t>(K),
        std::vector<uint8_t>(row_len, 0));

    for (int i = 0; i < K; i++) {
        for (size_t j = 0; j < row_len; j++) {
            uint8_t acc = 0;
            for (int k = 0; k < K_prime; k++) {
                uint8_t coef = V[static_cast<size_t>(i)][static_cast<size_t>(k)];
                uint8_t datum = data[static_cast<size_t>(k)][j];
                acc ^= field.mul(coef, datum);
            }
            encoded[static_cast<size_t>(i)][j] = acc;
        }
    }
    return encoded;
}

std::vector<std::vector<uint8_t>> ErasureCode::gaussian_solve(
    const GF256& field,
    std::vector<std::vector<uint8_t>> A,
    std::vector<std::vector<uint8_t>> B,
    int n, int m) {

    // Forward elimination with partial pivoting
    for (int col = 0; col < n; col++) {
        // Find pivot
        int pivot = -1;
        for (int row = col; row < n; row++) {
            if (A[static_cast<size_t>(row)][static_cast<size_t>(col)] != 0) {
                pivot = row;
                break;
            }
        }
        if (pivot == -1) {
            // Singular
            return {};
        }
        if (pivot != col) {
            std::swap(A[static_cast<size_t>(pivot)], A[static_cast<size_t>(col)]);
            std::swap(B[static_cast<size_t>(pivot)], B[static_cast<size_t>(col)]);
        }

        uint8_t pivot_val = A[static_cast<size_t>(col)][static_cast<size_t>(col)];
        uint8_t pivot_inv = field.inverse(pivot_val);

        // Normalize pivot row
        for (int c = col; c < n; c++) {
            A[static_cast<size_t>(col)][static_cast<size_t>(c)] =
                field.mul(A[static_cast<size_t>(col)][static_cast<size_t>(c)], pivot_inv);
        }
        for (int c = 0; c < m; c++) {
            B[static_cast<size_t>(col)][static_cast<size_t>(c)] =
                field.mul(B[static_cast<size_t>(col)][static_cast<size_t>(c)], pivot_inv);
        }

        // Eliminate other rows
        for (int row = 0; row < n; row++) {
            if (row == col) continue;
            uint8_t factor = A[static_cast<size_t>(row)][static_cast<size_t>(col)];
            if (factor == 0) continue;
            for (int c = col; c < n; c++) {
                A[static_cast<size_t>(row)][static_cast<size_t>(c)] ^=
                    field.mul(factor, A[static_cast<size_t>(col)][static_cast<size_t>(c)]);
            }
            for (int c = 0; c < m; c++) {
                B[static_cast<size_t>(row)][static_cast<size_t>(c)] ^=
                    field.mul(factor, B[static_cast<size_t>(col)][static_cast<size_t>(c)]);
            }
        }
    }
    return B; // B is now the solution
}

std::vector<std::vector<uint8_t>> ErasureCode::decode(
    int K, int K_prime,
    const std::vector<int>& surviving_indices,
    const std::vector<std::vector<uint8_t>>& surviving_data) {

    if (K <= 0 || K_prime <= 0 || K < K_prime) {
        throw std::invalid_argument("ErasureCode::decode: need 1 <= K' <= K");
    }
    if (K_prime > max_dimension() || K > max_dimension()) {
        throw std::invalid_argument(
            "ErasureCode::decode: K and K' must be <= 255 for GF(2^8)");
    }
    if (static_cast<int>(surviving_indices.size()) != K_prime) {
        throw std::invalid_argument(
            "ErasureCode::decode: need exactly K' surviving indices");
    }
    if (static_cast<int>(surviving_data.size()) != K_prime) {
        throw std::invalid_argument(
            "ErasureCode::decode: need exactly K' surviving data rows");
    }

    // Validate indices
    for (int i = 0; i < K_prime; i++) {
        int idx = surviving_indices[static_cast<size_t>(i)];
        if (idx < 0 || idx >= K) {
            throw std::invalid_argument(
                "ErasureCode::decode: surviving index out of range");
        }
    }
    // Check for duplicates
    auto sorted_idx = surviving_indices;
    std::sort(sorted_idx.begin(), sorted_idx.end());
    for (int i = 1; i < K_prime; i++) {
        if (sorted_idx[static_cast<size_t>(i)] == sorted_idx[static_cast<size_t>(i - 1)]) {
            throw std::invalid_argument(
                "ErasureCode::decode: duplicate surviving indices");
        }
    }

    // Determine row length
    size_t row_len = 0;
    if (!surviving_data.empty()) {
        row_len = surviving_data[0].size();
        for (const auto& row : surviving_data) {
            if (row.size() != row_len) {
                throw std::invalid_argument(
                    "ErasureCode::decode: all surviving rows must have same length");
            }
        }
    }

    GF256 field;

    // Build the full Vandermonde matrix (K x K')
    auto V_full = Vandermonde::build(field, K, K_prime);

    // Extract the K' x K' submatrix corresponding to surviving rows
    // A_sub[i][j] = V_full[surviving_indices[i]][j]
    std::vector<std::vector<uint8_t>> A_sub(
        static_cast<size_t>(K_prime),
        std::vector<uint8_t>(static_cast<size_t>(K_prime), 0));
    for (int i = 0; i < K_prime; i++) {
        int row_idx = surviving_indices[static_cast<size_t>(i)];
        for (int j = 0; j < K_prime; j++) {
            A_sub[static_cast<size_t>(i)][static_cast<size_t>(j)] =
                V_full[static_cast<size_t>(row_idx)][static_cast<size_t>(j)];
        }
    }

    // B = surviving_data (K' x row_len), copied
    std::vector<std::vector<uint8_t>> B(
        static_cast<size_t>(K_prime),
        std::vector<uint8_t>(row_len, 0));
    for (int i = 0; i < K_prime; i++) {
        for (size_t j = 0; j < row_len; j++) {
            B[static_cast<size_t>(i)][j] =
                surviving_data[static_cast<size_t>(i)][j];
        }
    }

    // Solve A_sub * x = B for x (the original data, K' x row_len)
    auto solution = gaussian_solve(field, A_sub, B, K_prime,
                                   static_cast<int>(row_len));
    if (solution.empty()) {
        throw std::runtime_error(
            "ErasureCode::decode: system is singular; cannot decode");
    }
    return solution;
}

} // namespace qkdpp
