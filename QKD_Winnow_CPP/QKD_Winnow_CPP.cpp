#include <iostream>
#include <numeric>
#include <vector>
#include <string>
#include <random>
#include <chrono>
#include <armadillo>
#include "BS_thread_pool.hpp"
#include <indicators/cursor_control.hpp>
#include <indicators/progress_bar.hpp>

#define DEBUG_WINNOW false

using namespace std;
using namespace arma;
using namespace chrono;
using namespace indicators;


const int TRIAL_NUMBER = 100;
const int THREADS_NUMBER = 16;
const bool SHUFFLE_MODE = true;
const int SIFTED_KEY_LENGTH = 10240;
const int INITIAL_SYNDROME_POWER = 3;
const string RESULT_FILE_PATH = "C:\\Users\\AGENT\\Desktop\\QKD_Winnow_CPP\\results\\";
const vector<double> BER = { 0.01, 0.03, 0.05, 0.06, 0.08, 0.1, 0.12, 0.15 };
const vector<vector<int>> COMBINATION_ELEMENTS = 
{   { 0, 1, 2, 3, 4 },
    { 0, 1, 2, 3, 4 },
    { 0, 1, 2, 3, 4 },
    { 0, 1, 2, 3, 4 },
    { 0, 1, 2, 3, 4 }
};

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
    double error_probability{};
};

void _PRINT_ARRAY(int* bit_array, size_t ba_size, size_t block_size) {
    for (size_t i = 0; i < ba_size; i++)
    {
        if (i % block_size == 0 && i != 0)
        {
            cout << " ";
        }
        cout << bit_array[i];
        
    }
    cout << endl;
}

string get_trial_combination_string(vector<int> combination) {
    string comb_str = "(";
    for (size_t i = 0; i < combination.size(); i++)
    {
        comb_str += to_string(combination[i]);
        if (i < combination.size() - 1)
        {
            comb_str += ",";
        }
        else
        {
            comb_str += ")";
        }
    }
    return comb_str;
}

bool write_file(vector<test_result>& data, string path) {
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

vector<vector<int>> cartesian_product(vector<vector<int>>& combinations) {
    auto product = [](long long a, vector<int>& b) { return a * b.size(); };
    const long long combination_number = accumulate(combinations.begin(), combinations.end(), 1LL, product);
    vector<vector<int>> result(combination_number, vector<int>(combinations.size()));
    for (long long n = 0; n < combination_number; ++n) {
        lldiv_t q{ n, 0 };
        for (long long i = combinations.size() - 1; 0 <= i; --i) {
            q = div(q.quot, combinations[i].size());
            result[n][i] = combinations[i][q.rem];
        }
    }
    return result;
}

vector<test_combination> prepare_combinations(vector<vector<int>> trial_elements, vector<double> bit_error_rates) {
    vector<vector<int>> trial_combinations = cartesian_product(trial_elements);
    vector<test_combination> combinations(trial_combinations.size() * bit_error_rates.size());
    int test_number = 0;
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

void generate_random_bit_array(size_t length, int* output_random_bit_array) {
    // Initialize a Mersenne Twister engine for random number generation
    mt19937 generator(time(0));
    uniform_int_distribution<int> distribution(0, 1); // Distribution for 0 or 1

    // Generate random bits and fill the vector
    for (int i = 0; i < length; ++i) {
        output_random_bit_array[i] = distribution(generator);
    }
}

void introduce_errors(int* bit_array, size_t array_length, double error_probability, int* output_bit_array_with_errors) {
    mt19937 generator(time(0));
    uniform_real_distribution<double> distribution(0.0, 1.0);

    for (int i = 0; i < array_length; ++i) {

        output_bit_array_with_errors[i] = bit_array[i] ^ (distribution(generator) < error_probability);
    }
}

void discard_bits_for_parity_check(int*source_bit_array, size_t source_array_length, int* destination_bit_array, size_t syndrome_power) {
    int source_block_size = (int)pow(2, syndrome_power);
    int destination_block_size = source_block_size - 1;
    for (size_t i = 0, j = 0; i < source_array_length; i += source_block_size, j += destination_block_size) {
        copy(source_bit_array + i + 1, source_bit_array + i + source_block_size, destination_bit_array + j);
    }
}

void discard_bits_for_syndrome(int* source_bit_block, int* destination_bit_block, vector<int>& discarded_bit_positions) {
    int destination_current_start = 0;
    for (size_t i = 1; i < discarded_bit_positions.size() - 1; i++)
    {
        destination_current_start += discarded_bit_positions[i] - discarded_bit_positions[i - 1] - 1;
        copy(source_bit_block + discarded_bit_positions[i] + 1, source_bit_block + discarded_bit_positions[i+1], destination_bit_block + destination_current_start);
    }
}

bool block_parity_check(int* bit_block, size_t block_length) {
    return accumulate(bit_block, bit_block + block_length, 0) % 2 == 0;
}

void calculate_error_positions(int* alice_bit_array, int* bob_bit_array, size_t array_length, int* output_error_positions) {
    for (size_t i = 0; i < array_length; i++)
    {
        output_error_positions[i] = alice_bit_array[i] ^ bob_bit_array[i];
    }
}

void shuffle_array_bits(int* alice_bit_array, int* bob_bit_array, size_t array_length, int seed) {
    mt19937 rng(seed);
    shuffle(alice_bit_array, alice_bit_array + array_length, rng);
    rng.seed(seed);
    shuffle(bob_bit_array, bob_bit_array + array_length, rng);
}

umat calculate_Hamming_hash_function(int syndrome_power) {
    int block_length = (int)pow(2, syndrome_power) - 1;

    umat hash_matrix(syndrome_power, block_length);
    for (int i = 0; i < hash_matrix.n_rows; i++)
    {
        for (int j = 0; j < hash_matrix.n_cols; j++)
        {
            hash_matrix(i, j) = (int)floor((j + 1) / pow(2, i)) % 2;
        }
    }
    return hash_matrix;
}

void calculate_syndrome(int* bit_block, size_t block_length, umat& hash_matrix, int* output_syndrome) {
    ucolvec column_vector(block_length);
    for (size_t i = 0; i < block_length; i++)
    {
        column_vector(i) = bit_block[i];
    }
    ucolvec syndrome = hash_matrix * column_vector;
 
    for (size_t i = 0; i < syndrome.size(); i++)
    {
        output_syndrome[i] = syndrome(i) % 2;
    }
}

void correct_error(int* bit_block, int* first_syndrome, int* second_syndrome, size_t syndrome_power) {
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

size_t winnow(int* alice_bit_array, int* bob_bit_array, size_t array_length, size_t syndrome_power, int* output_alice_bit_array, int* output_bob_bit_array) {
    
    size_t block_len = (int)pow(2, syndrome_power);
    int blocks_cnt = array_length / block_len;
    umat hash_mat = calculate_Hamming_hash_function(syndrome_power);

    vector<int> diff_par_blocks;
    diff_par_blocks.reserve(blocks_cnt);
    for (size_t i = 0; i < array_length; i += block_len)
    {
        if (block_parity_check(alice_bit_array + i, block_len) != block_parity_check(bob_bit_array + i, block_len))
        {
            diff_par_blocks.push_back((int)(i / block_len));
        }
    }
    #if DEBUG_WINNOW
    for (size_t i = 0; i < diff_par_blocks.size(); i++)
    {
        cout << diff_par_blocks[i] << ' ';
    }
    cout << endl;
    #endif

    size_t priv_amp_arr_len = array_length - (int)array_length / block_len;
    int* alice_priv_amp = new int[priv_amp_arr_len];
    int* bob_priv_amp = new int[priv_amp_arr_len];
    discard_bits_for_parity_check(alice_bit_array, array_length, alice_priv_amp, syndrome_power);
    discard_bits_for_parity_check(bob_bit_array, array_length, bob_priv_amp, syndrome_power);

    #if DEBUG_WINNOW
    _PRINT_ARRAY(alice_bit_array, array_length, block_len);
    _PRINT_ARRAY(bob_bit_array, array_length, block_len);
    cout << endl;
    _PRINT_ARRAY(alice_priv_amp, priv_amp_arr_len, block_len-1);
    _PRINT_ARRAY(bob_priv_amp, priv_amp_arr_len, block_len-1);
    cout << endl;
    #endif

    block_len -= 1;
    int* alice_syn = new int[syndrome_power];
    int* bob_syn = new int[syndrome_power];
    for (size_t i = 0; i < diff_par_blocks.size(); i++)
    {
        calculate_syndrome(alice_priv_amp + (diff_par_blocks[i] * block_len), block_len, hash_mat, alice_syn);
        calculate_syndrome(bob_priv_amp + (diff_par_blocks[i] * block_len), block_len, hash_mat, bob_syn);
        correct_error(alice_priv_amp + (diff_par_blocks[i] * block_len), alice_syn, bob_syn, syndrome_power);
    }
    delete[] alice_syn;
    delete[] bob_syn;
    #if DEBUG_WINNOW
    _PRINT_ARRAY(alice_priv_amp, priv_amp_arr_len, block_len);
    _PRINT_ARRAY(bob_priv_amp, priv_amp_arr_len, block_len);
    cout << endl;
    #endif 

    size_t remain_bits_cnt = block_len - syndrome_power;
    size_t out_arr_len = priv_amp_arr_len - diff_par_blocks.size() * syndrome_power;

    vector<int> disc_bit_pos(syndrome_power + 1);
    for (size_t i = 0; i < disc_bit_pos.size(); i++)
    {
        disc_bit_pos[i] = (int)(pow(2, i) - 1);
    }

    size_t j = 0;
    size_t copy_delta = 0;
    size_t dest_cur_pos = 0;
    size_t diff_par_len = diff_par_blocks.size();
    for (size_t i = 0; i < priv_amp_arr_len; )
    {
        if (j < diff_par_len && (i / block_len == diff_par_blocks[j]))
        {
            discard_bits_for_syndrome(alice_priv_amp + i, output_alice_bit_array + dest_cur_pos, disc_bit_pos);
            discard_bits_for_syndrome(bob_priv_amp + i, output_bob_bit_array + dest_cur_pos, disc_bit_pos);
            dest_cur_pos += remain_bits_cnt;
            i += block_len;
            j++;
        }
        else if (diff_par_len == 0) 
        {
            copy(alice_priv_amp, alice_priv_amp + priv_amp_arr_len, output_alice_bit_array);
            copy(bob_priv_amp, bob_priv_amp + priv_amp_arr_len, output_bob_bit_array);
            i += priv_amp_arr_len;
        }
        else if (j >= diff_par_len)
        {
            copy_delta = priv_amp_arr_len - i;
            copy(alice_priv_amp + i, alice_priv_amp + i + copy_delta, output_alice_bit_array + dest_cur_pos);
            copy(bob_priv_amp + i, bob_priv_amp + i + copy_delta, output_bob_bit_array + dest_cur_pos);
            i += copy_delta;
        }
        else
        {
            copy_delta = diff_par_blocks[j] * block_len - i;
            copy(alice_priv_amp + i, alice_priv_amp + i + copy_delta, output_alice_bit_array + dest_cur_pos);
            copy(bob_priv_amp + i, bob_priv_amp + i + copy_delta, output_bob_bit_array + dest_cur_pos);
            dest_cur_pos += copy_delta;
            i += copy_delta;
        }
    }
    #if DEBUG_WINNOW
    _PRINT_ARRAY(output_alice_bit_array, out_arr_len, block_len);
    _PRINT_ARRAY(output_bob_bit_array, out_arr_len, block_len);
    cout << endl;
    #endif

    delete[] alice_priv_amp;
    delete[] bob_priv_amp;

    return out_arr_len;
}

size_t run_trial(int* alice_bit_array, int* bob_bit_array, size_t array_length, vector<int>& trial_combination, bool shuffle_bits, int* output_alice_bit_array, int* output_bob_bit_array){
    int seed = 0;
    size_t block_length;
    size_t discarded_bits_number = 0;
    size_t trimmed_array_length = array_length;
    size_t syndrome_power = INITIAL_SYNDROME_POWER;
    int* current_alice_bit_array = new int[array_length];
    int* current_bob_bit_array = new int[array_length];

    memcpy(current_alice_bit_array, alice_bit_array, array_length * sizeof(int));
    memcpy(current_bob_bit_array, bob_bit_array, array_length * sizeof(int));
    for (size_t i = 0; i < trial_combination.size(); i++)
    {
        block_length = (int)pow(2, syndrome_power);
        for (size_t j = 0; j < trial_combination[i]; j++)
        {
            discarded_bits_number = trimmed_array_length % block_length;
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

test_result run_test(test_combination combination) {
    int errors_number = 0;
    size_t output_array_length = 0;
    double mean_final_error = 0;
    double mean_remaining_fraction = 0;
    int* alice_bit_array = new int[SIFTED_KEY_LENGTH];
    int* bob_bit_array = new int[SIFTED_KEY_LENGTH];
    int* output_alice_bit_array = new int[SIFTED_KEY_LENGTH] {};
    int* output_bob_bit_array = new int[SIFTED_KEY_LENGTH] {};
    int* error_positions_array = new int[SIFTED_KEY_LENGTH];

    for (size_t i = 0; i < TRIAL_NUMBER; i++)
    {
        generate_random_bit_array(SIFTED_KEY_LENGTH, alice_bit_array);
        introduce_errors(alice_bit_array, SIFTED_KEY_LENGTH, combination.error_probability, bob_bit_array);
        output_array_length = run_trial(alice_bit_array, bob_bit_array, SIFTED_KEY_LENGTH, combination.trial_combination, SHUFFLE_MODE, output_alice_bit_array, output_bob_bit_array);
        
        calculate_error_positions(output_alice_bit_array, output_bob_bit_array, output_array_length, error_positions_array);
        errors_number = accumulate(error_positions_array, error_positions_array + output_array_length, 0);
        mean_final_error += (double)errors_number / (double)output_array_length;
        mean_remaining_fraction += (double)output_array_length / (double)SIFTED_KEY_LENGTH;
    }
    mean_final_error = mean_final_error / (double)TRIAL_NUMBER;
    mean_remaining_fraction = mean_remaining_fraction / (double)TRIAL_NUMBER;

    delete[] alice_bit_array;
    delete[] bob_bit_array;
    delete[] output_alice_bit_array;
    delete[] output_bob_bit_array;
    delete[] error_positions_array;

    test_result result;
    result.test_number = combination.test_number;
    result.trial_combination = combination.trial_combination;
    result.error_probability = combination.error_probability;
    result.mean_final_error = mean_final_error;
    result.mean_remaining_fraction = mean_remaining_fraction;
    return result;
}

vector<test_result> run_simulation(vector<test_combination>& combinations) {
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

    int iteration = 1;
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
    srand(time(0));

    vector<test_combination> combinations = prepare_combinations(COMBINATION_ELEMENTS, BER);
    vector<test_result> result = run_simulation(combinations);
    write_file(result, RESULT_FILE_PATH);
    cout << endl;
    //for (size_t i = 0; i < result.size(); i++)
    //{
    //    for (size_t j = 0; j < result[i].size(); j++)
    //    {
    //        cout << result[i][j] <<" ";
    //    }
    //    cout << endl;
    //}


    //vector<int> comb = { 4,4,4,4,4 };

    //auto start_time = steady_clock::now();
    //for (size_t i = 0; i < 100; i++)
    //{
    //    test_result tst_res = run_test(comb, 0.15, 0);
    //}
    //auto end_time = steady_clock::now();
    //duration<double> execution_time = end_time - start_time;
    //cout << "\n" << execution_time << endl; 
     return 0;
}
