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
    float QBER{};
};

struct test_result
{
    int test_number{};
    std::vector<size_t> trial_combination{};
    double QBER{};

    double final_qber_mean{};
    double final_qber_std_dev{};
    double final_qber_min{};
    double final_qber_max{};

    double final_fraction_mean{};
    double final_fraction_std_dev{};
    double final_fraction_min{};
    double final_fraction_max{};
};

std::vector<std::vector<size_t>> cartesian_product(std::vector<std::vector<size_t>> trial_elements);
std::vector<test_combination> prepare_combinations(const std::vector<std::vector<size_t>>& trial_combinations, std::vector<double> bit_error_rates);
size_t run_trial(const int *const alice_bit_array, const int *const bob_bit_array, size_t array_length,
                 const std::vector<size_t> &trial_combination, bool shuffle_bits, int *const output_alice_bit_array, int *const output_bob_bit_array);
test_result run_test(const test_combination combination, size_t seed);
std::vector<test_result> run_simulation(const std::vector<test_combination> &combinations);
std::string get_trial_combination_string(const std::vector<size_t> &combination);
void write_file(const std::vector<test_result> &data, fs::path directory);
