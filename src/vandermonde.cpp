#include "qkdpp/vandermonde.hpp"
#include <stdexcept>

namespace qkdpp {

std::vector<std::vector<uint8_t>> Vandermonde::build(
    const GF256& field, int rows, int cols) {

    if (rows <= 0 || cols <= 0) {
        return {};
    }
    // The multiplicative group of GF(2^8) has order 255, and antilog() reduces
    // its exponent modulo 255.  Any index reaching 255 therefore repeats row 0
    // (or column 0), which would make the matrix rank-deficient and destroy the
    // MDS property.  Both dimensions are capped at 255 for that reason.
    if (cols > GF256::ORDER) {
        throw std::invalid_argument(
            "Vandermonde: cols must be <= 255 for GF(2^8)");
    }
    if (rows > GF256::ORDER) {
        throw std::invalid_argument(
            "Vandermonde: rows must be <= 255 for GF(2^8)");
    }

    // Row i: [alpha^(0*i), alpha^(1*i), alpha^(2*i), ..., alpha^((cols-1)*i)]
    // = [1, alpha^i, alpha^(2i), ..., alpha^((cols-1)*i)]
    std::vector<std::vector<uint8_t>> matrix(
        static_cast<size_t>(rows),
        std::vector<uint8_t>(static_cast<size_t>(cols), 0));

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            // alpha^(i*j) — computed in 64 bits: the int product i*j can
            // overflow for large dimensions even though both stay <= 255.
            long long exp = static_cast<long long>(i) * static_cast<long long>(j);
            matrix[static_cast<size_t>(i)][static_cast<size_t>(j)] =
                field.antilog(static_cast<int>(exp % GF256::ORDER));
        }
    }
    return matrix;
}

bool Vandermonde::rows_linearly_independent(
    const GF256& field,
    const std::vector<std::vector<uint8_t>>& matrix) {

    if (matrix.empty()) return true;

    // Reject ragged input: using row 0's width to index every row would read
    // (and in the row-reduction copy, write) out of bounds.
    const std::size_t m_sz = matrix[0].size();
    for (const auto& row : matrix) {
        if (row.size() != m_sz) {
            throw std::invalid_argument(
                "Vandermonde::rows_linearly_independent: matrix must be rectangular");
        }
    }

    int n = static_cast<int>(matrix.size());
    int m = static_cast<int>(m_sz);
    if (m == 0) return false;   // any non-empty set of zero-width rows is dependent

    // Build a copy for row reduction
    std::vector<std::vector<uint8_t>> A = matrix;

    int rank = 0;
    for (int col = 0; col < m && rank < n; col++) {
        // Find pivot row at or below `rank`
        int pivot = -1;
        for (int row = rank; row < n; row++) {
            if (A[static_cast<size_t>(row)][static_cast<size_t>(col)] != 0) {
                pivot = row;
                break;
            }
        }
        if (pivot == -1) continue;

        // Swap pivot row to `rank`
        if (pivot != rank) {
            std::swap(A[static_cast<size_t>(pivot)],
                      A[static_cast<size_t>(rank)]);
        }

        uint8_t pivot_val = A[static_cast<size_t>(rank)][static_cast<size_t>(col)];
        uint8_t pivot_inv = field.inverse(pivot_val);

        // Eliminate all other rows
        for (int row = 0; row < n; row++) {
            if (row == rank) continue;
            uint8_t factor = A[static_cast<size_t>(row)][static_cast<size_t>(col)];
            if (factor == 0) continue;
            uint8_t mult = field.mul(factor, pivot_inv);
            for (int c = col; c < m; c++) {
                A[static_cast<size_t>(row)][static_cast<size_t>(c)] ^=
                    field.mul(mult, A[static_cast<size_t>(rank)][static_cast<size_t>(c)]);
            }
        }
        rank++;
    }

    return rank == n;
}

} // namespace qkdpp
