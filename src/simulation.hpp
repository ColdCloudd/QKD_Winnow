#pragma once
#include <random>
#include <vector>
#include <numeric>
#include <cstdlib>

#include <BS_thread_pool.hpp>
#include <indicators/progress_bar.hpp>
#include <indicators/cursor_control.hpp>

#include "config.hpp"
#include "bit_array_operations.hpp"
#include "winnow_algorithm.hpp"

struct test_combination
{
    int test_number{};
    std::vector<size_t> trial_combination{};
    float error_probability{};
};

struct test_result
{
    int test_number{};
    std::vector<size_t> trial_combination{};
    double error_probability{};
    double mean_final_error{};
    double mean_remaining_fraction{};
};

std::vector<std::vector<size_t>> cartesian_product(std::vector<std::vector<size_t>> trial_elements);
std::vector<test_combination> prepare_combinations(std::vector<std::vector<size_t>> trial_elements, std::vector<double> bit_error_rates);
size_t run_trial(const int *const alice_bit_array, const int *const bob_bit_array, size_t array_length,
                 const std::vector<size_t> &trial_combination, bool shuffle_bits, int *const output_alice_bit_array, int *const output_bob_bit_array);
test_result run_test(const test_combination combination, size_t seed);
std::vector<test_result> run_simulation(const std::vector<test_combination> &combinations);
std::string get_trial_combination_string(const std::vector<size_t> &combination);
void write_file(const std::vector<test_result> &data, fs::path directory);
