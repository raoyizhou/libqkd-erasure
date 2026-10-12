#ifndef QKDPP_ERASURE_CODE_HPP
#define QKDPP_ERASURE_CODE_HPP

#include <cstdint>
#include <vector>

namespace qkdpp {

class GF256; // forward declaration

// [K, K'] Maximum Distance Separable (MDS) erasure code over GF(2^8).
//
// Encoding:   data (K' rows of equal length)  --V-->  encoded (K rows)
//   where V is a K x K' Vandermonde matrix.
//
// Decoding:   given any K' surviving rows out of the K encoded rows,
//   solve the K' x K' system to recover the original data.
//
// This implements the erasure coding described in §V-I/§V-J of the
// MP-QKD-TP protocol paper.  With K' = K/2, up to half the rows can be
// lost and the data is still recoverable.
class ErasureCode {
public:
    // Encode K' data rows into K encoded rows.
    // Each data row must have the same length.
    // K must be >= K' >= 1.
    static std::vector<std::vector<uint8_t>> encode(
        int K, int K_prime,
        const std::vector<std::vector<uint8_t>>& data);

    // Decode (recover original data) from any K' surviving rows.
    // surviving_indices: which of the K original encoded rows survived
    //   (size must == K_prime, values in [0, K))
    // surviving_data:    the data of those surviving rows (each same length)
    // Returns: K' reconstructed data rows.
    static std::vector<std::vector<uint8_t>> decode(
        int K, int K_prime,
        const std::vector<int>& surviving_indices,
        const std::vector<std::vector<uint8_t>>& surviving_data);

private:
    // Gaussian elimination solver over GF(2^8): solves A * x = b
    // where A is n x n, b is n x m (columns of data).
    // Returns the solution x (n rows, each of length m) or empty on failure.
    static std::vector<std::vector<uint8_t>> gaussian_solve(
        const GF256& field,
        std::vector<std::vector<uint8_t>> A,         // n x n (copied, will be modified)
        std::vector<std::vector<uint8_t>> B,         // n x m (copied, will be modified)
        int n, int m);
};

} // namespace qkdpp

#endif // QKDPP_ERASURE_CODE_HPP
