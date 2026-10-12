cmake_minimum_required(VERSION 3.15)
project(libqkd-erasure VERSION 1.0.0 LANGUAGES CXX)

# --- Build settings -------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Architecture-appropriate optimization flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Enable PCLMUL if the compiler supports it (used for potential GF multiply
# acceleration; safe to omit on CPUs without the instruction).
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-mpclmul" COMPILER_SUPPORTS_PCLMUL)
if(COMPILER_SUPPORTS_PCLMUL)
    add_compile_options(-mpclmul)
endif()

# Warnings
add_compile_options(-Wall -Wextra -Wpedantic)

# --- Library target -------------------------------------------------------
add_library(qkdpp_erasure STATIC
    src/gf2_8.cpp
    src/vandermonde.cpp
    src/erasure_code.cpp
)
target_include_directories(qkdpp_erasure PUBLIC include)

# --- Test target ---------------------------------------------------------
enable_testing()

add_executable(qkdpp_tests tests/test_main.cpp)
target_link_libraries(qkdpp_tests PRIVATE qkdpp_erasure)

# Output test binary to build/bin/
set_target_properties(qkdpp_tests PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

add_test(NAME qkdpp_tests COMMAND qkdpp_tests)

# --- Install (optional) --------------------------------------------------
install(TARGETS qkdpp_erasure
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)
install(DIRECTORY include/ DESTINATION include)
