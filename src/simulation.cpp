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
            combinations[test_number].QBER = bit_error_rates[j];
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
    size_t last_incompl_block_length = 0;
    size_t padding_length = 0;
    size_t num_padding_bits_removed_by_winnow = 0;
    size_t result_array_length = array_length;
    size_t syndrome_length = CFG.INITIAL_SYNDROME_LENGTH;
    int *current_alice_bit_array = new int[array_length];
    int *current_bob_bit_array = new int[array_length];
    bool padding = false;
    winnow_result result{};

    memcpy(current_alice_bit_array, alice_bit_array, array_length * sizeof(int));
    memcpy(current_bob_bit_array, bob_bit_array, array_length * sizeof(int));
    for (size_t i = 0; i < trial_combination.size(); i++)
    {
        block_length = static_cast<size_t>(pow(2, syndrome_length));
        int **hash_mat = calculate_Hamming_hash_matrix(syndrome_length);
        for (size_t j = 0; j < trial_combination[i]; j++)
        {
            last_incompl_block_length = result_array_length % block_length; // Before each Winnow run, сalculates the length of the last block
            
            if(last_incompl_block_length > 0)   // Padding 
            {
                padding_length = block_length - last_incompl_block_length;
                if(result_array_length + padding_length >= array_length)    // If padding extends beyond the initial array, the last block is cut off
                {
                     result_array_length -= last_incompl_block_length;
                }
                else
                {
                    std::fill(current_alice_bit_array + result_array_length, current_alice_bit_array + (result_array_length + padding_length), 0);
                    std::fill(current_bob_bit_array + result_array_length, current_bob_bit_array + (result_array_length + padding_length), 0);
                    result_array_length += padding_length;
                    padding = true;
                }
            }
            
            result = winnow(current_alice_bit_array, current_bob_bit_array, result_array_length, syndrome_length, hash_mat, output_alice_bit_array, output_bob_bit_array);
            result_array_length = result.result_array_length;

            if(padding)
            {
                if(result.last_block_with_error)
                {
                    for (size_t k = 0; k < syndrome_length; k++)
                    {
                        if (static_cast<size_t>(pow(2, k)) >= last_incompl_block_length)
                        {
                            num_padding_bits_removed_by_winnow++;       // Number of padding bits that were removed by the winnow process 
                        }
                    } 
                }
                result_array_length -= (padding_length - num_padding_bits_removed_by_winnow);
                num_padding_bits_removed_by_winnow = 0;
                padding = false;
            }
            
            if (shuffle_bits)
            {
                shuffle_array_bits(output_alice_bit_array, output_bob_bit_array, result_array_length, seed);
                seed++;
            }
            memcpy(current_alice_bit_array, output_alice_bit_array, result_array_length * sizeof(int));
            memcpy(current_bob_bit_array, output_bob_bit_array, result_array_length * sizeof(int));
        }

        for (size_t i = 0; i < syndrome_length; ++i)
        {
            delete[] hash_mat[i];
        }
        delete[] hash_mat;
        syndrome_length++;
    }
    delete[] current_alice_bit_array;
    delete[] current_bob_bit_array;

    return result_array_length;
}

// Runs the experiment with the given TRIAL_NUMBER- times combination and calculates
// the final average key error rate and the final average key fraction
test_result run_test(const test_combination combination, size_t seed)
{
    size_t errors_number = 0;
    size_t output_array_length = 0;
    int *alice_bit_array = new int[CFG.SIFTED_KEY_LENGTH];
    int *bob_bit_array = new int[CFG.SIFTED_KEY_LENGTH];
    int *output_alice_bit_array = new int[CFG.SIFTED_KEY_LENGTH]{};
    int *output_bob_bit_array = new int[CFG.SIFTED_KEY_LENGTH]{};
    int *error_positions_array = new int[CFG.SIFTED_KEY_LENGTH];

    std::vector<double> final_qbers(CFG.TRIALS_NUMBER); 
    std::vector<double> final_fractions(CFG.TRIALS_NUMBER); 

    // Pseudo-random number generator
    std::mt19937 prng(seed);

    for (size_t i = 0; i < CFG.TRIALS_NUMBER; i++)
    {
        generate_random_bit_array(prng, CFG.SIFTED_KEY_LENGTH, alice_bit_array);
        introduce_errors(prng, alice_bit_array, CFG.SIFTED_KEY_LENGTH, combination.QBER, bob_bit_array);
        output_array_length = run_trial(alice_bit_array, bob_bit_array, CFG.SIFTED_KEY_LENGTH, combination.trial_combination, CFG.SHUFFLE_MODE, output_alice_bit_array, output_bob_bit_array);

        calculate_error_positions(output_alice_bit_array, output_bob_bit_array, output_array_length, error_positions_array);
        errors_number = std::accumulate(error_positions_array, error_positions_array + output_array_length, 0);
        final_qbers[i] = static_cast<double>(errors_number) / static_cast<double>(output_array_length);
        final_fractions[i] = static_cast<double>(output_array_length) / static_cast<double>(CFG.SIFTED_KEY_LENGTH);
    }

    delete[] alice_bit_array;
    delete[] bob_bit_array;
    delete[] output_alice_bit_array;
    delete[] output_bob_bit_array;
    delete[] error_positions_array;

    double final_qber_mean = std::accumulate(final_qbers.begin(), final_qbers.end(), 0.) / static_cast<double>(CFG.TRIALS_NUMBER);
    auto final_qber_minmax = std::minmax_element(final_qbers.begin(), final_qbers.end());
    double final_qber_min = *final_qber_minmax.first;
    double final_qber_max = *final_qber_minmax.second;
    double final_qber_variance = std::accumulate(final_qbers.begin(), final_qbers.end(), 0., 
        [final_qber_mean](double acc, double val) {
            return acc + std::pow(val - final_qber_mean, 2);
        });
    double final_qber_std_dev = sqrt(final_qber_variance / CFG.TRIALS_NUMBER);

    double final_fraction_mean = std::accumulate(final_fractions.begin(), final_fractions.end(), 0.) / static_cast<double>(CFG.TRIALS_NUMBER);
    auto final_fraction_minmax = std::minmax_element(final_fractions.begin(), final_fractions.end());
    double final_fraction_min = *final_fraction_minmax.first;
    double final_fraction_max = *final_fraction_minmax.second;
    double final_fraction_variance = std::accumulate(final_fractions.begin(), final_fractions.end(), 0., 
        [final_fraction_mean](double acc, double val) {
            return acc + std::pow(val - final_fraction_mean, 2);
        });
    double final_fraction_std_dev = sqrt(final_fraction_variance / CFG.TRIALS_NUMBER);

    test_result result;
    result.test_number = combination.test_number;
    result.trial_combination = combination.trial_combination;
    result.initial_qber = static_cast<double>(static_cast<size_t>(CFG.SIFTED_KEY_LENGTH * combination.QBER)) / CFG.SIFTED_KEY_LENGTH; // Exact QBER in the key.

    result.final_qber_mean = final_qber_mean;
    result.final_qber_std_dev = final_qber_std_dev;
    result.final_qber_min = final_qber_min;
    result.final_qber_max = final_qber_max;

    result.final_fraction_mean = final_fraction_mean;
    result.final_fraction_std_dev = final_fraction_std_dev;
    result.final_fraction_min = final_fraction_min;
    result.final_fraction_max = final_fraction_max;

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

// Returns the run represented by the sequence (example: 0;2;1 - 0 passes with block size 8, then 2 passes with block size 16 and 1 pass with block size 32).
std::string get_num_pass_with_block_size_sequence_string(const std::vector<size_t> &combination)
{
    std::string seq_str = "";
    for (size_t i = 0; i < combination.size(); i++)
    {
        seq_str += std::to_string(combination[i]);
        if (i < combination.size() - 1)
        {
            seq_str += ";";
        }
    }
    return seq_str;
}

// Returns a block size header (example: N=8;N=16;N=32 ...).
std::string get_header_block_size_string(const std::vector<size_t> &combination)
{
    std::string header_str = "";
    size_t syndrome_length = CFG.INITIAL_SYNDROME_LENGTH;
    for (size_t i = 0; i < combination.size(); i++)
    {
        header_str += "N=";
        header_str += std::to_string(static_cast<size_t>(pow(2, syndrome_length)));
        if (i < combination.size() - 1)
        {
            header_str += ";";
        }
        syndrome_length++;
    }
    return header_str;
}

// Records the results of the simulation in a ".csv" format file
void write_file(const std::vector<test_result> &data, fs::path directory)
{
    try
    {
        if (!fs::exists(directory))
        {
            fs::create_directories(directory);
        }
        std::string filename = "winnow(trial_num=" + std::to_string(CFG.TRIALS_NUMBER) + ",shuff_mode=" + std::to_string(CFG.SHUFFLE_MODE) 
        + ",seed=" + std::to_string(CFG.SIMULATION_SEED) + ",key_length=" + std::to_string(CFG.SIFTED_KEY_LENGTH) + ").csv";
        fs::path result_file_path = directory / filename;

        std::fstream fout;
        fout.open(result_file_path, std::ios::out | std::ios::trunc);
        fout << "№;TRIAL_COMBINATION;"<< get_header_block_size_string(data[0].trial_combination) << ";INITIAL_QBER;FINAL_QBER_MEAN;FINAL_QBER_STD_DEV;FINAL_QBER_MIN;FINAL_QBER_MAX;" 
             << "FINAL_FRACTION_MEAN;FINAL_FRACTION_STD_DEV;FINAL_FRACTION_MIN;FINAL_FRACTION_MAX\n";
        for (size_t i = 0; i < data.size(); i++)
        {
            fout << data[i].test_number << ";" << get_trial_combination_string(data[i].trial_combination) << ";"<< get_num_pass_with_block_size_sequence_string(data[i].trial_combination) << ";"
                 << data[i].initial_qber << ";" << data[i].final_qber_mean << ";" << data[i].final_qber_std_dev << ";" << data[i].final_qber_min << ";" << data[i].final_qber_max << ";"
                 << data[i].final_fraction_mean << ";" << data[i].final_fraction_std_dev << ";" << data[i].final_fraction_min << ";" << data[i].final_fraction_max << "\n";
        }
        fout.close();
    }
    catch (const std::exception &e)
    {
        fmt::print(stderr, fg(fmt::color::red),"An error occurred while writing to the file.\n");
        throw;
    }
}