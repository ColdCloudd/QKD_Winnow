#pragma once 
#include <random>
#include <algorithm>

#include <XoshiroCpp.hpp>

void generate_random_bit_array(XoshiroCpp::Xoshiro256PlusPlus &prng, size_t length, int *const output_random_bit_array);
void introduce_errors(XoshiroCpp::Xoshiro256PlusPlus &prng, const int *const bit_array, size_t array_length, float QBER,
                      int *const output_bit_array_with_errors);
void shuffle_array_bits(int *const alice_bit_array, int *const bob_bit_array, size_t array_length, int seed);
void calculate_error_positions(const int *const alice_bit_array, const int *const bob_bit_array, size_t array_length,
                               int *const output_error_positions);                      