#ifndef QKDPP_VANDERMONDE_HPP
#define QKDPP_VANDERMONDE_HPP

#include <cstdint>
#include <vector>
#include "gf2_8.hpp"

namespace qkdpp {

// Constructs a (rows x cols) Vandermonde matrix over GF(2^8).
// Row i: [1, alpha^i, alpha^(2i), ..., alpha^((cols-1)*i)]
// where alpha is the primitive element of the GF256 field.
//
// This matrix is used as the MDS encoding matrix:
//   encoded = V * data
// where data has K' rows and V is K x K'.
class Vandermonde {
public:
    // Build a (rows x cols) Vandermonde matrix.
    // Requires cols <= 256 for the GF(2^8) field to support enough
    // distinct powers (any row count is fine).
    static std::vector<std::vector<uint8_t>> build(
        const GF256& field, int rows, int cols);

    // Verify that the first `rank` rows are linearly independent
    // (used in testing; also validates MDS property for square submatrices).
    static bool rows_linearly_independent(
        const GF256& field,
        const std::vector<std::vector<uint8_t>>& matrix);
};

} // namespace qkdpp

#endif // QKDPP_VANDERMONDE_HPP
