#include "simulation.hpp"

// Computes combinations of runs, based on given elements as the Cartesian product of vectors
std::vector<std::vector<size_t>> cartesian_product(std::vector<std::vector<size_t>> trial_elements)
{
    auto product = [](long long a, std::vector<size_t> &b)
    { return a * b.size(); };
    const long long combination_number = accumulate(trial_elements.begin(), trial_elements.end(), 1LL, product);
    std::vector<std::vector<size_t>> result(combination_number, std::vector<size_t>(trial_elements.size()));
    for (long long n = 0; n < combination_number; ++n)
    {
        std::lldiv_t q{n, 0};
        for (long long i = trial_elements.size() - 1; 0 <= i; --i)
        {
            q = std::div(q.quot, trial_elements[i].size());
            result[n][i] = trial_elements[i][q.rem];
        }
    }
    return result;
}

// Generates combinations of runs consisting of the run number, combinations that determine
// the number of Winnow runs with a given block size, and the QBER probability
std::vector<test_combination> prepare_combinations(const std::vector<std::vector<size_t>>& trial_combinations, std::vector<double> bit_error_rates)
{
    std::vector<test_combination> combinations(trial_combinations.size() * bit_error_rates.size());
    size_t test_number = 0;
    for (size_t i = 0; i < trial_combinations.size(); i++)
    {
        for (size_t j = 0; j < bit_error_rates.size(); j++)
        {
            combinations[test_number].test_number = test_number;
            combinations[test_number].trial_combination = trial_combinations[i];
            combinations[test_number].error_probability = bit_error_rates[j];
            test_number++;
        }
    }
    return combinations;
}

// Runs the Winnow algorithm sequentially several times with different block sizes
size_t run_trial(const int *const alice_bit_array, const int *const bob_bit_array, size_t array_length,
                 const std::vector<size_t> &trial_combination, bool shuffle_bits, int *const output_alice_bit_array, int *const output_bob_bit_array)
{
    size_t seed = CFG.SIMULATION_SEED;
    size_t block_length = 0;
    size_t discarded_bits_number = 0;
    size_t trimmed_array_length = array_length;
    size_t syndrome_power = CFG.INITIAL_SYNDROME_POWER;
    int *current_alice_bit_array = new int[array_length];
    int *current_bob_bit_array = new int[array_length];

    memcpy(current_alice_bit_array, alice_bit_array, array_length * sizeof(int));
    memcpy(current_bob_bit_array, bob_bit_array, array_length * sizeof(int));
    for (size_t i = 0; i < trial_combination.size(); i++)
    {
        block_length = static_cast<size_t>(pow(2, syndrome_power));
        for (size_t j = 0; j < trial_combination[i]; j++)
        {
            discarded_bits_number = trimmed_array_length % block_length; // Before each Winnow run, trims bit arrays to be a multiple of the current block length
            trimmed_array_length -= discarded_bits_number;
            trimmed_array_length = winnow(current_alice_bit_array, current_bob_bit_array, trimmed_array_length, syndrome_power, output_alice_bit_array, output_bob_bit_array);
            if (shuffle_bits)
            {
                shuffle_array_bits(output_alice_bit_array, output_bob_bit_array, trimmed_array_length, seed);
                seed++;
            }
            memcpy(current_alice_bit_array, output_alice_bit_array, trimmed_array_length * sizeof(int));
            memcpy(current_bob_bit_array, output_bob_bit_array, trimmed_array_length * sizeof(int));
        }
        syndrome_power++;
    }
    delete[] current_alice_bit_array;
    delete[] current_bob_bit_array;

    return trimmed_array_length;
}

// Runs the experiment with the given TRIAL_NUMBER- times combination and calculates
// the final average key error rate and the final average key fraction
test_result run_test(const test_combination combination, size_t seed)
{
    size_t errors_number = 0;
    size_t output_array_length = 0;
    double mean_final_error = 0;
    double mean_remaining_fraction = 0;
    int *alice_bit_array = new int[CFG.SIFTED_KEY_LENGTH];
    int *bob_bit_array = new int[CFG.SIFTED_KEY_LENGTH];
    int *output_alice_bit_array = new int[CFG.SIFTED_KEY_LENGTH]{};
    int *output_bob_bit_array = new int[CFG.SIFTED_KEY_LENGTH]{};
    int *error_positions_array = new int[CFG.SIFTED_KEY_LENGTH];

    // Pseudo-random number generator
    std::mt19937 prng(seed);

    for (size_t i = 0; i < CFG.TRIALS_NUMBER; i++)
    {
        generate_random_bit_array(prng, CFG.SIFTED_KEY_LENGTH, alice_bit_array);
        introduce_errors(prng, alice_bit_array, CFG.SIFTED_KEY_LENGTH, combination.error_probability, bob_bit_array);
        output_array_length = run_trial(alice_bit_array, bob_bit_array, CFG.SIFTED_KEY_LENGTH, combination.trial_combination, CFG.SHUFFLE_MODE, output_alice_bit_array, output_bob_bit_array);

        calculate_error_positions(output_alice_bit_array, output_bob_bit_array, output_array_length, error_positions_array);
        errors_number = std::accumulate(error_positions_array, error_positions_array + output_array_length, 0);
        mean_final_error += static_cast<double>(errors_number) / static_cast<double>(output_array_length);
        mean_remaining_fraction += static_cast<double>(output_array_length) / static_cast<double>(CFG.SIFTED_KEY_LENGTH);
    }
    mean_final_error = mean_final_error / static_cast<double>(CFG.TRIALS_NUMBER);
    mean_remaining_fraction = mean_remaining_fraction / static_cast<double>(CFG.TRIALS_NUMBER);

    delete[] alice_bit_array;
    delete[] bob_bit_array;
    delete[] output_alice_bit_array;
    delete[] output_bob_bit_array;
    delete[] error_positions_array;

    test_result result;
    result.test_number = combination.test_number;
    result.trial_combination = combination.trial_combination;
    result.error_probability = static_cast<double>(combination.error_probability);
    result.mean_final_error = mean_final_error;
    result.mean_remaining_fraction = mean_remaining_fraction;
    return result;
}

// Distributes all combinations of the experiment evenly across the CPU threads and runs it
std::vector<test_result> run_simulation(const std::vector<test_combination> &combinations)
{
    using namespace indicators;

    std::vector<test_result> results(combinations.size());
    BS::thread_pool pool(CFG.THREADS_NUMBER);

    indicators::show_console_cursor(false);
    indicators::ProgressBar bar{
        option::BarWidth{50},
        option::Start{" ["},
        option::Fill{"="},
        option::Lead{">"},
        option::Remainder{"-"},
        option::End{"]"},
        option::PrefixText{"PROGRESS"},
        option::ForegroundColor{Color::green},
        option::ShowElapsedTime{true},
        option::ShowRemainingTime{true},
        option::FontStyles{std::vector<FontStyle>{FontStyle::bold}},
        option::MaxProgress{combinations.size()}};

    size_t iteration = 0;
    std::mt19937 prng(CFG.SIMULATION_SEED);
    std::uniform_int_distribution<size_t> distribution(0, std::numeric_limits<size_t>::max());
    pool.detach_loop<size_t>(0, combinations.size(),
                             [&combinations, &results, &prng, &distribution, &bar, &iteration](size_t i)
                             {
                                 bar.set_option(option::PostfixText{
                                     std::to_string(iteration) + "/" + std::to_string(combinations.size())});
                                 bar.tick();
                                 iteration++;

                                 results[i] = run_test(combinations[i], distribution(prng));
                             });
    pool.wait();
    indicators::show_console_cursor(true);

    return results;
}

// Returns the combination as a python tuple in string format
std::string get_trial_combination_string(const std::vector<size_t> &combination)
{
    std::string comb_str = "(";
    for (size_t i = 0; i < combination.size(); i++)
    {
        comb_str += std::to_string(combination[i]);
        if (i < combination.size() - 1)
        {
            comb_str += ", ";
        }
    }
    comb_str += ")";
    return comb_str;
}

// Records the results of the simulation in a ".csv" format file
void write_file(const std::vector<test_result> &data, fs::path directory)
{
    try
    {
        std::string filename = "winnow(trial_num=" + std::to_string(CFG.TRIALS_NUMBER) + ",shuff_mode=" + std::to_string(CFG.SHUFFLE_MODE) + ",seed=" + std::to_string(CFG.SIMULATION_SEED) + ").csv";
        fs::path result_file_path = directory / filename;

        std::fstream fout;
        fout.open(result_file_path, std::ios::out | std::ios::trunc);
        fout << "№;TRIAL_COMBINATION;INITIAL_QBER;MEAN_FINAL_QBER;MEAN_FINAL_FRACTION\n";
        for (size_t i = 0; i < data.size(); i++)
        {
            fout << data[i].test_number << ";" << get_trial_combination_string(data[i].trial_combination) << ";" << data[i].error_probability << ";"
                 << data[i].mean_final_error << ";" << data[i].mean_remaining_fraction << "\n";
        }
        fout.close();
    }
    catch (const std::exception &e)
    {
        fmt::print(stderr, fg(fmt::color::red),"An error occurred while writing to the file.\n");
        throw;
    }
}