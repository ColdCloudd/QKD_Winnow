#pragma once
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>

#include <fmt/core.h>
#include <fmt/color.h>
#include <fmt/ranges.h>

#include "config.hpp"
#include "utils.hpp"

void discard_bits_for_parity_check(const int *const source_bit_array, const size_t &source_array_length,
                                   int *const destination_bit_array, const size_t &syndrome_power);
void discard_bits_for_syndrome(const int *const source_bit_block, int *const destination_bit_block,
                               const std::vector<int> &discarded_bit_positions);
bool calculate_block_parity(const int *const bit_block, const size_t &block_length);
int **calculate_Hamming_hash_matrix(size_t syndrome_power);
void calculate_syndrome(const int *const bit_block, const size_t &syndrome_power, const size_t &block_length,
                        const int *const *hash_matrix, int *const output_syndrome);                               
void correct_error(int *const bit_block, const int *const first_syndrome, const int *const second_syndrome,
                   const size_t &syndrome_power);
size_t winnow(int *const alice_bit_array, int *const bob_bit_array, size_t array_length, size_t syndrome_power,
              int *const output_alice_bit_array, int *const output_bob_bit_array);                   

