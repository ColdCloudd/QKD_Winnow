#include <iostream>
#include <numeric>
#include <vector>
#include <string>
#include <random>
#include <chrono>
#include <armadillo>

using namespace std;
using namespace chrono;
using namespace arma;

int INITIAL_SYNDROME_POWER = 3;
int SIFTED_KEY_LENGTH = 10240;
vector<float> BER = { 0.01, 0.03, 0.05, 0.06, 0.08, 0.1, 0.12, 0.15 };
vector<vector<int>> COMBINATION_ELEMENTS {{ 0, 1, 2 },
                                          { 0, 1, 2 },
                                          { 0, 1, 2 }};
int TRIAL_NUMBER = 10;
int PROCESSES_NUMBER = 12;
bool SHUFFLE_MODE = true;
string RESULT_FILE_PATH = "C:\\Users\\AGENT\\Desktop\\QKD_Winnow_CPP\\results\\winnow_results.csv";

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

void generate_random_bit_array(size_t length, int* output_random_bit_array) {
    // Initialize a Mersenne Twister engine for random number generation
    mt19937 generator(time(0));
    uniform_int_distribution<int> distribution(0, 1); // Distribution for 0 or 1

    // Generate random bits and fill the vector
    for (int i = 0; i < length; ++i) {
        output_random_bit_array[i] = distribution(generator);
    }
}

void introduce_errors(int* bit_array, size_t array_length, float initial_error_probability, int* output_bit_array_with_errors) {
    mt19937 generator(time(0));
    uniform_real_distribution<double> distribution(0.0, 1.0);

    for (int i = 0; i < array_length; ++i) {

        output_bit_array_with_errors[i] = bit_array[i] ^ (distribution(generator) < initial_error_probability);
    }
}

void discard_bits_for_parity_check(int*source_bit_array, size_t source_array_length, int* destination_bit_array, size_t syndrome_power) {
    int source_block_size = pow(2, syndrome_power);
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

void shuffle_bits(int* alice_bit_array, int* bob_bit_array, size_t array_length, int seed) {
    mt19937 rng(seed);
    shuffle(alice_bit_array, alice_bit_array + array_length, rng);
    rng.seed(seed);
    shuffle(bob_bit_array, bob_bit_array + array_length, rng);
}

umat calculate_Hamming_hash_function(int syndrome_power) {
    int block_length = pow(2, syndrome_power) - 1;

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
        error_bit_position += (first_syndrome[i] ^ second_syndrome[i]) * pow(2, i);
    }
    if (error_bit_position >= 0)
    {
        bit_block[error_bit_position] = !bit_block[error_bit_position];
    }
}

int winnow(int* alice_bit_array, int* bob_bit_array, size_t array_length, size_t syndrome_power, int* output_alice_bit_array, int* output_bob_bit_array) {
    size_t block_length = pow(2, syndrome_power);
    int blocks_count = array_length / block_length;
    umat hash_matrix = calculate_Hamming_hash_function(syndrome_power);

    vector<int> b_odd_parity_block_numbers;
    b_odd_parity_block_numbers.reserve(blocks_count);
    for (size_t i = 0; i < array_length; i += block_length)
    {
        if (block_parity_check(alice_bit_array + i, block_length) != block_parity_check(bob_bit_array + i, block_length))
        {
            b_odd_parity_block_numbers.push_back(i / block_length);
        }
    }
    /*for (size_t i = 0; i < b_odd_parity_block_numbers.size(); i++)
    {
        cout << b_odd_parity_block_numbers[i] << ' ';
    }
    cout << endl;*/

    size_t amp1_array_size = array_length - (int)array_length / block_length;
    int* alice_privacy_amp_1 = new int[amp1_array_size];
    int* bob_privacy_amp_1 = new int[amp1_array_size];
    discard_bits_for_parity_check(alice_bit_array, array_length, alice_privacy_amp_1, syndrome_power);
    discard_bits_for_parity_check(bob_bit_array, array_length, bob_privacy_amp_1, syndrome_power);
    /*_PRINT_ARRAY(alice_privacy_amp_1, amp1_array_size, 7);
    _PRINT_ARRAY(bob_privacy_amp_1, amp1_array_size, 7);*/

    block_length -= 1;
    int* alice_syndrome = new int[syndrome_power];
    int* bob_syndrome = new int[syndrome_power];
    for (size_t i = 0; i < b_odd_parity_block_numbers.size(); i++)
    {
        calculate_syndrome(alice_privacy_amp_1 + (b_odd_parity_block_numbers[i] * block_length), block_length, hash_matrix, alice_syndrome);
        calculate_syndrome(bob_privacy_amp_1 + (b_odd_parity_block_numbers[i] * block_length), block_length, hash_matrix, bob_syndrome);
        correct_error(alice_privacy_amp_1 + (b_odd_parity_block_numbers[i] * block_length), alice_syndrome, bob_syndrome, syndrome_power);
    }
    delete[] alice_syndrome;
    delete[] bob_syndrome;
    /*_PRINT_ARRAY(alice_privacy_amp_1, amp1_array_size, 7);
    _PRINT_ARRAY(bob_privacy_amp_1, amp1_array_size, 7);*/
    

    int remaining_bits_count = block_length - syndrome_power;
    size_t output_array_size = amp1_array_size - b_odd_parity_block_numbers.size() * syndrome_power;

    vector<int> discarded_bit_positions(syndrome_power + 1);
    for (size_t i = 0; i < discarded_bit_positions.size(); i++)
    {
        discarded_bit_positions[i] = pow(2, i) - 1;
    }

    int j = 0;
    int copy_delta = 0;
    int dest_cur_pos = 0;
    int b_odd_size = b_odd_parity_block_numbers.size();
    for (size_t i = 0; i < amp1_array_size; )
    {
        if (j < b_odd_size && (i / block_length == b_odd_parity_block_numbers[j]))
        {
            discard_bits_for_syndrome(alice_privacy_amp_1 + i, output_alice_bit_array + dest_cur_pos, discarded_bit_positions);
            discard_bits_for_syndrome(bob_privacy_amp_1 + i, output_bob_bit_array + dest_cur_pos, discarded_bit_positions);
            dest_cur_pos += remaining_bits_count;
            i += block_length;
            j++;
        }
        else if (b_odd_size == 0) 
        {
            copy(alice_privacy_amp_1, alice_privacy_amp_1 + amp1_array_size, output_alice_bit_array);
            copy(bob_privacy_amp_1, bob_privacy_amp_1 + amp1_array_size, output_bob_bit_array);
            i += amp1_array_size;
        }
        else if (j >= b_odd_size)
        {
            copy_delta = amp1_array_size - i;
            copy(alice_privacy_amp_1 + i, alice_privacy_amp_1 + i + copy_delta, output_alice_bit_array + dest_cur_pos);
            copy(bob_privacy_amp_1 + i, bob_privacy_amp_1 + i + copy_delta, output_bob_bit_array + dest_cur_pos);
            i += copy_delta;
        }
        else
        {
            copy_delta = b_odd_parity_block_numbers[j] * block_length - i;
            copy(alice_privacy_amp_1 + i, alice_privacy_amp_1 + i + copy_delta, output_alice_bit_array + dest_cur_pos);
            copy(bob_privacy_amp_1 + i, bob_privacy_amp_1 + i + copy_delta, output_bob_bit_array + dest_cur_pos);
            dest_cur_pos += copy_delta;
            i += copy_delta;
        }
    }
   
    /*_PRINT_ARRAY(alice_privacy_amp_1, amp1_array_size, 7);
    _PRINT_ARRAY(alice_privacy_amp_2, amp2_array_size, 7);*/
    delete[] alice_privacy_amp_1;
    delete[] bob_privacy_amp_1;

    return output_array_size;
}

int main() {
    srand(time(0));
    
    /*vector<vector<int>> result = cartesian_product(COMBINATION_ELEMENTS);
    for (size_t i = 0; i < result.size(); i++)
    {
        for (size_t j = 0; j < result[i].size(); j++)
        {
            cout << result[i][j] <<" ";
        }
        cout << endl;
    }*/
    int* random_bit_array = new int[SIFTED_KEY_LENGTH];
    int* bit_array_with_errors = new int[SIFTED_KEY_LENGTH];
    
    //_PRINT_ARRAY(random_bit_array, SIFTED_KEY_LENGTH, 8);
    //_PRINT_ARRAY(bit_array_with_errors, SIFTED_KEY_LENGTH, 8);
    int* output_alice_bit_array = new int[SIFTED_KEY_LENGTH];
    int* output_bob_bit_array = new int[SIFTED_KEY_LENGTH];

    auto start_time = steady_clock::now();
    for (size_t i = 0; i < 10000; i++)
    {
        generate_random_bit_array(SIFTED_KEY_LENGTH, random_bit_array);
        introduce_errors(random_bit_array, SIFTED_KEY_LENGTH, 0.15, bit_array_with_errors);
        winnow(random_bit_array, bit_array_with_errors, SIFTED_KEY_LENGTH, INITIAL_SYNDROME_POWER, output_alice_bit_array, output_bob_bit_array);
    }
    delete[] output_alice_bit_array;
    delete[] output_bob_bit_array;
    auto end_time = steady_clock::now();
    duration<double> execution_time = end_time - start_time;
    cout << "\n" << execution_time << endl; 

    //int dest_cur_pos = 0;
    //int b_odd_size = b_odd_parity_block_numbers.size();
    //int j = 0;
    //for (size_t i = 0; i < amp1_array_size; i += block_length)
    //{
    //    if (j < b_odd_size && (i / block_length == b_odd_parity_block_numbers[j]))
    //    {
    //        discard_bits_for_syndrome(alice_privacy_amp_1 + i, block_length, alice_privacy_amp_2 + dest_cur_pos, discarded_bit_positions);
    //        dest_cur_pos += remaining_bits_count;
    //        j++;
    //    }
    //    else
    //    {
    //        copy(alice_privacy_amp_1 + i, alice_privacy_amp_1 + i + block_length, alice_privacy_amp_2 + dest_cur_pos);
    //        dest_cur_pos += block_length;
    //    }
    //}
    
    
     return 0;
}
