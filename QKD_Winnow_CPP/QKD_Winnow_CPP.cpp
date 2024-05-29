#include <vector>
#include <string>
#include <random>
#include <chrono>
#include <numeric>
#include <iostream>
#include <fstream>

#include "BS_thread_pool.hpp"
#include <indicators/progress_bar.hpp>
#include <indicators/cursor_control.hpp>

#define DEBUG_MAIN false
#define DEBUG_WINNOW false

using namespace std;
using namespace chrono;
using namespace indicators;

//------------------------CONFIG PARAMETERS------------------------
// Number of runs with one combination
const int TRIAL_NUMBER = 100;

// Number of threads for parallelizing runs
const int THREADS_NUMBER = 16;

// Bit shuffling between protocol iterations 
const bool SHUFFLE_MODE = true;

// Initial key size
const int SIFTED_KEY_LENGTH = 10240;

// Initial power of syndrome. Determines the length of the block (2^3, then 2^4, etc.)
const int INITIAL_SYNDROME_POWER = 3;

// Average error rate in the key
const vector<double> QBER = { 0.01, 0.03, 0.05, 0.06, 0.08, 0.1, 0.12, 0.15 };

// Folder path for saving ".csv" files with the results of the experiment 
const string RESULT_FILE_PATH = "C:\\Users\\AGENT\\Desktop\\QKD_Winnow_CPP\\results\\";

// The numbers in the first vector determine the number of Winnow runs with a block length of 2^3, 
// the numbers in the second vector determine the number of runs with a block length of 2^4, and so on.
// This is necessary to compose all combinations of interest by calculating the Cartesian product.
const vector<vector<int>> COMBINATION_ELEMENTS = 
{   { 0, 1, 2, 3, 4 },
    { 0, 1, 2, 3, 4 },
    { 0, 1, 2, 3, 4 },
    { 0, 1, 2, 3, 4 },
    { 0, 1, 2, 3, 4 }
};
//------------------------CONFIG PARAMETERS------------------------

struct test_result
{
    int test_number{};
    vector<int> trial_combination{};
    double error_probability{};
    double mean_final_error{};
    double mean_remaining_fraction{};
};

struct test_combination
{
    int test_number{};
    vector<int> trial_combination{};
    float error_probability{};
};

void dbg_print_array(const int* const bit_array, size_t array_length, size_t block_length) {
    for (size_t i = 0; i < array_length; i++)
    {
        if (i % block_length == 0 && i != 0)
        {
            std::cout << " ";
        }
        std::cout << bit_array[i];
    }
    std::cout << endl;
}

// Returns the combination as a python tuple in string format
string get_trial_combination_string(const vector<int>& combination) {
    string comb_str = "(";
    for (size_t i = 0; i < combination.size(); i++)
    {
        comb_str += to_string(combination[i]);
        if (i < combination.size() - 1)
        {
            comb_str += ", ";
        }
        else
        {
            comb_str += ")";
        }
    }
    return comb_str;
}

// Records the results of the experiments in a ".csv" format file
bool write_file(const vector<test_result>& data, string path) {
    try
    {
        string filename = path + "winnow_res_cpp" + (string)"(trial_num=" + to_string(TRIAL_NUMBER)
            + (string)",shuff_mode=" + to_string(SHUFFLE_MODE) + ")" + (string)".csv";
        fstream fout;
        fout.open(filename, ios::out | ios::trunc);
        for (size_t i = 0; i < data.size(); i++)
        {
            fout << data[i].test_number << ";" << get_trial_combination_string(data[i].trial_combination) << ";" << data[i].error_probability << ";"
                << data[i].mean_final_error << ";" << data[i].mean_remaining_fraction << "\n";
        }
        fout.close();
        return true;
    }
    catch (const std::exception& ex)
    {
        printf("Error occured (write_file): %s", ex.what());
        return false;
    }
}

// Computes combinations of runs, based on given elements as the Cartesian product of vectors
vector<vector<int>> cartesian_product(vector<vector<int>> trial_elements) {
    auto product = [](long long a, vector<int>& b) { return a * b.size(); };
    const long long combination_number = accumulate(trial_elements.begin(), trial_elements.end(), 1LL, product);
    vector<vector<int>> result(combination_number, vector<int>(trial_elements.size()));
    for (long long n = 0; n < combination_number; ++n) {
        lldiv_t q{ n, 0 };
        for (long long i = trial_elements.size() - 1; 0 <= i; --i) {
            q = div(q.quot, trial_elements[i].size());
            result[n][i] = trial_elements[i][q.rem];
        }
    }
    return result;
}

// Generates combinations of runs consisting of the run number, combinations that determine 
// the number of Winnow runs with a given block size, and the QBER probability
vector<test_combination> prepare_combinations(vector<vector<int>> trial_elements, vector<double> bit_error_rates) {
    vector<vector<int>> trial_combinations = cartesian_product(trial_elements);
    vector<test_combination> combinations(trial_combinations.size() * bit_error_rates.size());
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

// Generates Alice's key
void generate_random_bit_array(mt19937& prng, size_t length, int* const output_random_bit_array) {
    uniform_int_distribution<int> distribution(0, 1);

    // Generate random bits and fill the vector
    for (int i = 0; i < length; ++i) {
        output_random_bit_array[i] = distribution(prng);
    }
}

// Generates Bob's key by making errors in Alice's key with a given QBER probability (Uniform distribution)
void introduce_errors(mt19937& prng, const int* const bit_array, size_t array_length, float error_probability, int* const output_bit_array_with_errors) {
    size_t num_errors = static_cast<size_t>(array_length * error_probability);
    if (num_errors == 0)
    {
        copy(bit_array, bit_array + array_length, output_bit_array_with_errors);
    }
    else
    {
        size_t* error_positions = new size_t[array_length];
        for (size_t i = 0; i < array_length; ++i)
        {
            error_positions[i] = i;
        }

        shuffle(error_positions, error_positions + array_length, prng);
        copy(bit_array, bit_array + array_length, output_bit_array_with_errors);

        for (size_t i = 0; i < num_errors; ++i)
        {
            output_bit_array_with_errors[error_positions[i]] ^= 1;
        }

        delete[] error_positions;
    }
}

// Discards bits at the first position in the block for privacy amplification
void discard_bits_for_parity_check(const int* const source_bit_array, const size_t& source_array_length, int* const destination_bit_array, const size_t& syndrome_power) {
    size_t source_block_size = static_cast<size_t>(pow(2, syndrome_power));
    size_t destination_block_size = source_block_size - 1;
    for (size_t i = 0, j = 0; i < source_array_length; i += source_block_size, j += destination_block_size) {
        copy(source_bit_array + i + 1, source_bit_array + i + source_block_size, destination_bit_array + j);
    }
}

// Discards bits at positions 2^n, where n=(0,1, ... syndrome_power - 1) in the block for privacy amplification
void discard_bits_for_syndrome(const int* const source_bit_block, int* const destination_bit_block, const vector<int>& discarded_bit_positions) {
    int destination_current_start = 0;
    // Copying bits that are between the positions described in discarded_bit_positions 
    for (size_t i = 1; i < discarded_bit_positions.size() - 1; i++)
    {
        destination_current_start += discarded_bit_positions[i] - discarded_bit_positions[i - 1] - 1;
        copy(source_bit_block + discarded_bit_positions[i] + 1, source_bit_block + discarded_bit_positions[i+1], destination_bit_block + destination_current_start);
    }
}

// Calculates the parity bit for a block
bool calculate_block_parity(const int* const bit_block, const size_t& block_length) {
    return accumulate(bit_block, bit_block + block_length, 0) % 2 == 0;
}

// Calculates an array that consists of 0 and 1, where 1 denotes that Alice's and Bob's bits are different
void calculate_error_positions(const int* const alice_bit_array, const int* const bob_bit_array, size_t array_length, int* const output_error_positions) {
    for (size_t i = 0; i < array_length; i++)
    {
        output_error_positions[i] = alice_bit_array[i] ^ bob_bit_array[i];
    }
}

// Shuffles Alice's and Bob's bit arrays by seed
void shuffle_array_bits(int* const alice_bit_array, int* const bob_bit_array, size_t array_length, int seed) {
    mt19937 rng(seed);
    shuffle(alice_bit_array, alice_bit_array + array_length, rng);
    rng.seed(seed);
    shuffle(bob_bit_array, bob_bit_array + array_length, rng);
}

// Computes a matrix based on the Hamming hash function
int** calculate_Hamming_hash_matrix(size_t syndrome_power) {
    size_t block_length = static_cast<size_t>(pow(2, syndrome_power)) - 1;

    int** hash_matrix = new int* [syndrome_power];
    for (size_t i = 0; i < syndrome_power; ++i)
    {
        hash_matrix[i] = new int[block_length];
        for (size_t j = 0; j < block_length; ++j)
        {
            hash_matrix[i][j] = static_cast<int>(floor((j + 1) / pow(2, i))) % 2;
        }
    }
    return hash_matrix;
}

// Calculates the syndrome by multiplying the Hamming matrix by the bit block column vector
void calculate_syndrome(const int* const bit_block, const size_t& syndrome_power, const size_t& block_length, const int* const* hash_matrix, int* const output_syndrome) {
    int xor_sum;
    for (size_t i = 0; i < syndrome_power; i++)
    {
        xor_sum = 0;
        for (size_t j = 0; j < block_length; j++)
        {
            if (bit_block[j])
            {
                xor_sum ^= hash_matrix[i][j];
            }
        }
        output_syndrome[i] = xor_sum;
    }
}

// Calculates the error position in a block based on Alice's and Bob's syndromes and inverts this bit in Alice's block
void correct_error(int* const bit_block, const int* const first_syndrome, const int* const second_syndrome, const size_t& syndrome_power) {
    int error_bit_position = -1;
    for (size_t i = 0; i < syndrome_power; i++)
    {
        error_bit_position += (first_syndrome[i] ^ second_syndrome[i]) * (int)pow(2, i);
    }
    if (error_bit_position >= 0)
    {
        bit_block[error_bit_position] = !bit_block[error_bit_position];
    }
}

size_t winnow(int* const alice_bit_array, int* const bob_bit_array, size_t array_length, size_t syndrome_power, int* const output_alice_bit_array, int* const output_bob_bit_array) {
    size_t block_len = static_cast<size_t>(pow(2, syndrome_power));
    size_t blocks_cnt = array_length / block_len;
    int** hash_mat = calculate_Hamming_hash_matrix(syndrome_power);

    vector<int> diff_par_blocks;    // Contains the numbers of blocks whose parity bits did not match for Alice and Bob
    diff_par_blocks.reserve(blocks_cnt);
    for (size_t i = 0; i < array_length; i += block_len)
    {
        if (calculate_block_parity(alice_bit_array + i, block_len) != calculate_block_parity(bob_bit_array + i, block_len))
        {
            diff_par_blocks.push_back(static_cast<int>(i / block_len));
        }
    }

    #if DEBUG_WINNOW
    std::cout << "______________DEBUG_WINNOW______________" << endl;
    for (size_t i = 0; i < diff_par_blocks.size(); i++)
    {
        std::cout << diff_par_blocks[i] << ' ';
    }
    std::cout << endl;
    #endif

    size_t priv_amp_arr_len = array_length - static_cast<size_t>(array_length / block_len);
    int* alice_priv_amp = new int[priv_amp_arr_len];
    int* bob_priv_amp = new int[priv_amp_arr_len];
    // Privacy amplification by discarding the first bit in each block.
    discard_bits_for_parity_check(alice_bit_array, array_length, alice_priv_amp, syndrome_power);
    discard_bits_for_parity_check(bob_bit_array, array_length, bob_priv_amp, syndrome_power);

    #if DEBUG_WINNOW
    dbg_print_array(alice_bit_array, array_length, block_len);
    dbg_print_array(bob_bit_array, array_length, block_len);
    std::cout << endl;
    dbg_print_array(alice_priv_amp, priv_amp_arr_len, block_len-1);
    dbg_print_array(bob_priv_amp, priv_amp_arr_len, block_len-1);
    std::cout << endl;
    #endif

    block_len -= 1;
    // Alice and Bob syndromes
    int* alice_syn = new int[syndrome_power];
    int* bob_syn = new int[syndrome_power];
    for (size_t i = 0; i < diff_par_blocks.size(); i++)
    {
        // Calculation of syndromes for blocks with non-matching parity bits, followed by error correction
        calculate_syndrome(alice_priv_amp + (diff_par_blocks[i] * block_len), syndrome_power, block_len, hash_mat, alice_syn);
        calculate_syndrome(bob_priv_amp + (diff_par_blocks[i] * block_len), syndrome_power, block_len, hash_mat, bob_syn);
        correct_error(alice_priv_amp + (diff_par_blocks[i] * block_len), alice_syn, bob_syn, syndrome_power);
    }

    for (size_t i = 0; i < syndrome_power; ++i) {
        delete[] hash_mat[i];
    }
    delete[] hash_mat;
    delete[] alice_syn;
    delete[] bob_syn;
    

    #if DEBUG_WINNOW
    dbg_print_array(alice_priv_amp, priv_amp_arr_len, block_len);
    dbg_print_array(bob_priv_amp, priv_amp_arr_len, block_len);
    std::cout << endl;
    #endif 

    size_t remain_bits_cnt = block_len - syndrome_power;   // Number of remaining bits in blocks for which syndromes were calculated
    size_t out_arr_len = priv_amp_arr_len - diff_par_blocks.size() * syndrome_power;

    // Contains bounds that specify valid bits [from 2^0, 2^1, ... , 2^(syndrome_power-1)], 
    // and the last element (2^syndrome_power) which is the right boundary
    vector<int> disc_bit_pos(syndrome_power + 1);   
    for (size_t i = 0; i < disc_bit_pos.size(); i++)
    {
        disc_bit_pos[i] = (int)(pow(2, i) - 1);
    }

    size_t j = 0;               // Used to move through diff_par_blocks
    size_t copy_delta = 0;      // Specifies the number of bits that can be copied 
    size_t dest_cur_pos = 0;    // The current position in the destination array from which new bits can be inserted
    size_t diff_par_len = diff_par_blocks.size();
    for (size_t i = 0; i < priv_amp_arr_len; )
    {
        if (j < diff_par_len && (i / block_len == diff_par_blocks[j]))      // In a block whose number is in diff_par_blocks, bits are discarded
        {
            discard_bits_for_syndrome(alice_priv_amp + i, output_alice_bit_array + dest_cur_pos, disc_bit_pos);
            discard_bits_for_syndrome(bob_priv_amp + i, output_bob_bit_array + dest_cur_pos, disc_bit_pos);
            dest_cur_pos += remain_bits_cnt;
            i += block_len;
            j++;
        }
        else if (diff_par_len == 0)     // No errors in arrays -> copy the array completely
        {
            copy(alice_priv_amp, alice_priv_amp + priv_amp_arr_len, output_alice_bit_array);
            copy(bob_priv_amp, bob_priv_amp + priv_amp_arr_len, output_bob_bit_array);
            i += priv_amp_arr_len;
        }
        else if (j >= diff_par_len)     // Remaining blocks with no errors are copied
        {
            copy_delta = priv_amp_arr_len - i;
            copy(alice_priv_amp + i, alice_priv_amp + i + copy_delta, output_alice_bit_array + dest_cur_pos);
            copy(bob_priv_amp + i, bob_priv_amp + i + copy_delta, output_bob_bit_array + dest_cur_pos);
            i += copy_delta;
        }
        else
        {
            copy_delta = diff_par_blocks[j] * block_len - i;    // Blocks between two erroneous blocks are copied
            copy(alice_priv_amp + i, alice_priv_amp + i + copy_delta, output_alice_bit_array + dest_cur_pos);
            copy(bob_priv_amp + i, bob_priv_amp + i + copy_delta, output_bob_bit_array + dest_cur_pos);
            dest_cur_pos += copy_delta;
            i += copy_delta;
        }
    }
    #if DEBUG_WINNOW
    dbg_print_array(output_alice_bit_array, out_arr_len, block_len);
    dbg_print_array(output_bob_bit_array, out_arr_len, block_len);
    std::cout << "______________DEBUG_WINNOW______________\n" << endl;
    #endif

    delete[] alice_priv_amp;
    delete[] bob_priv_amp;

    return out_arr_len;
}

// Runs the Winnow algorithm sequentially several times with different block sizes
size_t run_trial(const int* const alice_bit_array, const int* const bob_bit_array, size_t array_length, const vector<int>& trial_combination, bool shuffle_bits, int* const output_alice_bit_array, int* const output_bob_bit_array){
    int seed = 0;
    size_t block_length = 0;
    size_t discarded_bits_number = 0;
    size_t trimmed_array_length = array_length;
    size_t syndrome_power = INITIAL_SYNDROME_POWER;
    int* current_alice_bit_array = new int[array_length];
    int* current_bob_bit_array = new int[array_length];

    memcpy(current_alice_bit_array, alice_bit_array, array_length * sizeof(int));
    memcpy(current_bob_bit_array, bob_bit_array, array_length * sizeof(int));
    for (size_t i = 0; i < trial_combination.size(); i++)
    {
        block_length = static_cast<size_t>(pow(2, syndrome_power));
        for (size_t j = 0; j < trial_combination[i]; j++)
        {
            discarded_bits_number = trimmed_array_length % block_length;    // Before each Winnow run, trims bit arrays to be a multiple of the current block length
            trimmed_array_length -= discarded_bits_number;
            trimmed_array_length = winnow(current_alice_bit_array, current_bob_bit_array, trimmed_array_length, syndrome_power, output_alice_bit_array, output_bob_bit_array);
            if (shuffle_bits) {
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
test_result run_test(const test_combination combination)  {
    size_t errors_number = 0;
    size_t output_array_length = 0;
    double mean_final_error = 0;
    double mean_remaining_fraction = 0;
    int* alice_bit_array = new int[SIFTED_KEY_LENGTH];
    int* bob_bit_array = new int[SIFTED_KEY_LENGTH];
    int* output_alice_bit_array = new int[SIFTED_KEY_LENGTH] {};
    int* output_bob_bit_array = new int[SIFTED_KEY_LENGTH] {};
    int* error_positions_array = new int[SIFTED_KEY_LENGTH];

    // Pseudo-random number generators
    mt19937 prng_1(time(nullptr));
    mt19937 prng_2(time(nullptr) + 666);

    for (size_t i = 0; i < TRIAL_NUMBER; i++)
    {
        generate_random_bit_array(prng_1, SIFTED_KEY_LENGTH, alice_bit_array);
        introduce_errors(prng_2, alice_bit_array, SIFTED_KEY_LENGTH, combination.error_probability, bob_bit_array);
        output_array_length = run_trial(alice_bit_array, bob_bit_array, SIFTED_KEY_LENGTH, combination.trial_combination, SHUFFLE_MODE, output_alice_bit_array, output_bob_bit_array);
        
        calculate_error_positions(output_alice_bit_array, output_bob_bit_array, output_array_length, error_positions_array);
        errors_number = accumulate(error_positions_array, error_positions_array + output_array_length, 0);
        mean_final_error += static_cast<double>(errors_number) / static_cast<double>(output_array_length);
        mean_remaining_fraction += static_cast<double>(output_array_length) / static_cast<double>(SIFTED_KEY_LENGTH);
    }
    mean_final_error = mean_final_error / static_cast<double>(TRIAL_NUMBER);
    mean_remaining_fraction = mean_remaining_fraction / static_cast<double>(TRIAL_NUMBER);

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
vector<test_result> run_all_experiments(const vector<test_combination>& combinations) {
    vector<test_result> results(combinations.size());
    BS::thread_pool pool(THREADS_NUMBER);

    show_console_cursor(false);
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
      option::MaxProgress{combinations.size()}
    };

    size_t iteration = 1;
    pool.detach_loop<size_t>(0, combinations.size(),
        [&combinations, &results, &bar, &iteration](size_t i)
        {
            results[i] = run_test(combinations[i]);
            bar.set_option(option::PostfixText{
                to_string(iteration) + "/" + to_string(combinations.size())
                });
            bar.tick();
            iteration++;
        });
    pool.wait();
    show_console_cursor(true);

    return results;
}

int main() {
    #if not DEBUG_MAIN
    srand(time(nullptr));

    vector<test_combination> combinations = prepare_combinations(COMBINATION_ELEMENTS, QBER);
    vector<test_result> result = run_all_experiments(combinations);
    std::cout << "All tests were completed successfully! \nThe results will be written to the directory: " << RESULT_FILE_PATH << endl;
    if (write_file(result, RESULT_FILE_PATH))
    {
        printf("The results were written to the file successfully!");
    }
    else 
    {
        printf("An error occurred while writing to the file");
    }

    #else
    int* a_array = new int[32] 
        {0, 1, 0, 1, 1, 1, 0, 0,
         1, 1, 0, 1, 0, 1, 1, 1,
         0, 1, 0, 1, 1, 1, 0, 0,
         0, 1, 1, 1, 0, 1, 0, 1};
    int* b_array = new int[32] 
       {0, 1, 0, 1, 1, 1, 0, 0,
        1, 1, 0, 1, 0, 1, 1, 1,
        0, 1, 0, 1, 1, 1, 0, 0,
        0, 1, 1, 1, 0, 1, 0, 1};

    b_array[10] = 1;    // block #2

    b_array[17] = 0;    // block #3
    b_array[23] = 1;

    b_array[24] = 1;    // block #4
    b_array[25] = 0;
    b_array[26] = 0;

    std::cout << "______________DEBUG_MAIN______________" << endl;
    dbg_print_array(a_array, 32, 8);
    dbg_print_array(b_array, 32, 8);
    std::cout << "_ _ _ _ _ _ _ DEBUG_MAIN _ _ _ _ _ _ _\n" << endl;

    int* a_out_array = new int[32];
    int* b_out_array = new int[32];
    size_t out_len = winnow(a_array, b_array, 32, INITIAL_SYNDROME_POWER, a_out_array, b_out_array);

    std::cout << "_ _ _ _ _ _ _ DEBUG_MAIN _ _ _ _ _ _ _" << endl;
    dbg_print_array(a_out_array, out_len, 100);
    dbg_print_array(b_out_array, out_len, 100);
    std::cout << "______________DEBUG_MAIN______________" << endl;
    #endif
     return 0;
}
