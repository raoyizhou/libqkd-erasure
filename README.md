# libqkd-erasure

GF(2^8) erasure code engine for the **MP-QKD-TP** (Multi-Path Quantum Key Distribution Transport Protocol).

Implements the `[K, K']` MDS (Maximum Distance Separable) erasure coding described in §V-I / §V-J of the protocol paper, using a Vandermonde matrix over GF(2^8) (the AES field, polynomial `0x11B`).

## Features

- **GF(2^8) field arithmetic** with log/antilog (Zech) tables for O(1) multiply/inverse
- **Vandermonde matrix** construction: `V[i][j] = alpha^(i*j)`
- **[K, K'] MDS erasure encode/decode**: tolerates up to `K - K'` row losses
- **Gaussian elimination** solver over GF(2^8) for arbitrary surviving row subsets
- With the paper's default `K = 4, K' = 2` deployment, **up to half the paths can be lost** and the data is still fully recoverable

## Project structure

```
libqkd-erasure/
  include/qkdpp/
    gf2_8.hpp          GF(2^8) field arithmetic (log/antilog tables)
    vandermonde.hpp     Vandermonde matrix construction + independence check
    erasure_code.hpp    [K, K'] MDS encode/decode interface
  src/
    gf2_8.cpp
    vandermonde.cpp
    erasure_code.cpp
  tests/
    test_main.cpp       Unit tests (GF arithmetic, Vandermonde, erasure encode/decode)
  CMakeLists.txt
  README.md
  .gitignore
```

## Building

### Prerequisites

- C++17 compiler (GCC 13+ recommended, e.g. MinGW-w64 GCC 13.2.0)
- CMake 3.15+ (optional; see the manual command below for a single-step build)

### Manual build (no CMake)

```sh
g++ -std=c++17 -O3 -mpclmul -Wall -Wextra -Iinclude \
    src/*.cpp tests/test_main.cpp \
    -o build/bin/qkdpp_tests.exe
```

Then run:

```sh
./build/bin/qkdpp_tests.exe
```

### CMake build

```sh
cmake -B build -S .
cmake --build build --config Release
ctest --test-dir build --output-on-failure
```

## API

### GF(2^8) — `qkdpp::GF256`

```cpp
GF256 field;
uint8_t sum = field.add(a, b);      // XOR
uint8_t prod = field.mul(a, b);     // O(1) via log/antilog tables
uint8_t inv = field.inverse(a);     // multiplicative inverse
uint8_t q = field.div(a, b);       // a / b
uint8_t e = field.pow(a, exponent); // a^e
```

### Erasure code — `qkdpp::ErasureCode`

```cpp
// Encode K' data rows into K encoded rows (K >= K')
std::vector<std::vector<uint8_t>> encode(
    int K, int K_prime,
    const std::vector<std::vector<uint8_t>>& data);

// Decode from any K' surviving rows
std::vector<std::vector<uint8_t>> decode(
    int K, int K_prime,
    const std::vector<int>& surviving_indices,
    const std::vector<std::vector<uint8_t>>& surviving_data);
```

### Usage example

```cpp
#include "qkdpp/erasure_code.hpp"
#include "qkdpp/gf2_8.hpp"

using namespace qkdpp;

// K = 4 paths, K' = 2 minimum shares (tolerates 2 losses)
std::vector<std::vector<uint8_t>> data = {
    {0x42, 0x01, 0xFF, 0xAA},
    {0x80, 0x7F, 0x00, 0x55}
};
auto encoded = ErasureCode::encode(4, 2, data);

// Simulate losing rows 0 and 1, keeping rows 2 and 3
auto recovered = ErasureCode::decode(4, 2, {2, 3},
    {encoded[2], encoded[3]});
// recovered == data
```

## Technical details

### Field

- Polynomial: `0x11B` (x^8 + x^4 + x^3 + x + 1) — the AES irreducible polynomial
- Primitive element: `alpha = 0x03`
- Precomputed `log`/`antilog` tables (255 entries each) provide O(1) multiply, inverse, divide, and power

### Vandermonde matrix

Row `i` of the `K x K'` encoding matrix is:

```
[ 1, alpha^i, alpha^(2i), ..., alpha^((K'-1)*i) ]
```

This is an MDS (Maximum Distance Separable) matrix: any `K' x K'` submatrix is non-singular, guaranteeing that any `K'` of the `K` encoded rows suffice to recover the original data.

### Erasure code

- **Encoding**: `encoded = V * data`, where `V` is `K x K'` and `data` is `K' x L`
- **Decoding**: given any `K'` surviving rows, extract the corresponding `K' x K'` submatrix of `V`, then solve `V_sub * x = surviving_data` via Gaussian elimination over GF(2^8)

### Performance notes

- The `-mpclmul` flag is enabled where supported for potential PCLMULQDQ acceleration of GF multiply (the current table-based implementation does not require it, but it is harmless and future-profs vectorized paths)
- Log/antilog tables give O(1) field operations; encoding is `O(K * K' * L)` and decoding is `O(K'^3 + K'^2 * L)`

## Testing

The test suite (`tests/test_main.cpp`) covers:

- **GF(2^8)**: multiply (table vs Russian-peasant), zero, identity, commutativity, associativity, distributivity, inverse (a * a^-1 = 1), double inverse, division, power, log/antilog round-trip
- **Vandermonde**: shape, first row/column all-ones, linear independence, any submatrix non-singular (MDS property)
- **Erasure encode**: manual hand-computed verification for `K=4, K'=2`
- **Erasure decode**:
  - No loss (first `K'` rows)
  - 1 row lost
  - `K - K'` rows lost (maximum tolerance)
  - All `C(K, K')` row subsets
  - 50 random trials with random data and random surviving subsets
  - Larger `K=8, K'=4` configuration matching the paper's "half rows lost" scenario
  - Multi-byte (256-byte) and edge cases (1-byte, 0-byte rows)
