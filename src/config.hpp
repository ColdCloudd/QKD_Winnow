#pragma once
#include <vector>
#include <fstream>
#include <filesystem>

#include <fmt/core.h>
#include <fmt/color.h>
#include <fmt/ranges.h>
#include <nlohmann/json.hpp>

using json = nlohmann::json;
namespace fs = std::filesystem;

struct config_data
{
    // Number of threads for parallelizing runs
    size_t THREADS_NUMBER{};

    // Number of runs with one combination
    size_t TRIALS_NUMBER{};

    // Seed of simulation
    size_t SIMULATION_SEED{};

    // Bit shuffling between winnow iterations
    bool SHUFFLE_MODE{};

    // Initial key size
    size_t SIFTED_KEY_LENGTH{};

    // Initial power of syndrome. Determines the length of the block (2^3, then 2^4, etc.)
    size_t INITIAL_SYNDROME_POWER{};

    // Average initial error rate in the Bob's key
    std::vector<double> QBER{};

    // The numbers in the first vector determine the number of Winnow runs with a block length of 2^3,
    // the numbers in the second vector determine the number of runs with a block length of 2^4, and so on.
    // This is necessary to compose all combinations of interest by calculating the Cartesian product.
    std::vector<std::vector<size_t>> COMBINATION_ELEMENTS{};

    // Using combinations defined in the configuration, instead of generating them from the COMBINATION_ELEMENTS set
    bool USE_SPECIFIED_COMBINATIONS{};

    // Combinations that determine the number of runs with a given block size, starting from 2^3
    std::vector<std::vector<size_t>> COMBINATIONS{};

    // Flag for controlling the output of debugging information during the execution of winnow algorithm
    bool TRACE_WINNOW{};
};

extern config_data CFG;

config_data get_config_data(fs::path config_path);