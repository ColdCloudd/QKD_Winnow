#include <cstdlib>
#include <iostream>
#include <numeric>
#include <vector>
#include <string>
#include <ctime>
#include <random>
#include <chrono>
#include <armadillo>
#include <valarray>

using namespace std;
using namespace chrono;
using namespace arma;

int INITIAL_SYNDROME_POWER = 3;
int SIFTED_KEY_LENGTH = 24;
vector<float> BER = { 0.01, 0.03, 0.05, 0.06, 0.08, 0.1, 0.12, 0.15 };
vector<vector<int>> COMBINATION_ELEMENTS {{ 0, 1, 2 },
                                                { 0, 1, 2 },
                                                { 0, 1, 2 }};
int TRIAL_NUMBER = 10;
int PROCESSES_NUMBER = 12;
bool SHUFFLE_MODE = true;
string RESULT_FILE_PATH = "C:\\Users\\AGENT\\Desktop\\QKD_Winnow_NumPy\\results\\winnow_results.csv";

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

int* generate_random_bit_array(size_t length) {
    // Initialize a Mersenne Twister engine for random number generation
    mt19937 generator(time(0));
    uniform_int_distribution<int> distribution(0, 1); // Distribution for 0 or 1

    int* random_bit_array = new int[length];

    // Generate random bits and fill the vector
    for (int i = 0; i < length; ++i) {
        random_bit_array[i] = distribution(generator);
    }

    return random_bit_array;
}

int* introduce_errors(int* bit_array, size_t ba_size, float initial_error_probability) {
    mt19937 generator(time(0));
    uniform_real_distribution<double> distribution(0.0, 1.0);

    int* bit_array_with_errors = new int[ba_size];
    for (int i = 0; i < ba_size; ++i) {

        bit_array_with_errors[i] = bit_array[i] ^ (distribution(generator) < initial_error_probability);
    }

    return bit_array_with_errors;
}

void discard_bits_for_parity_check(int*source_bit_array, size_t sba_size, int*destination_bit_array, size_t dba_size, size_t syndrome_power) {
    int source_block_size = pow(2, syndrome_power);
    int destination_block_size = source_block_size - 1;
    int j = 0;
    for (size_t i = 0, j = 0; i < sba_size && j < dba_size; i += source_block_size, j += destination_block_size) {
        copy(source_bit_array + i + 1, source_bit_array + i + source_block_size, destination_bit_array + j);
    }
}

void discard_bits_for_syndrome(int* source_bit_block, size_t& sbb_size, int* destination_bit_block, vector<int>& discarded_bit_positions) {
    int destination_current_start = 0;
    for (size_t i = 1; i < discarded_bit_positions.size() - 1; i++)
    {
        destination_current_start += discarded_bit_positions[i] - discarded_bit_positions[i - 1] - 1;
        copy(source_bit_block + discarded_bit_positions[i] + 1, source_bit_block + discarded_bit_positions[i+1], destination_bit_block + destination_current_start);
    }
}

bool block_parity_check(int* bit_block, int block_size) {
    return accumulate(bit_block, bit_block + block_size, 0) % 2 == 0;
}

int* calculate_error_positions(int* alice_bit_array, int* bob_bit_array, size_t ba_size) {
    int* error_positions = new int[ba_size];
    for (size_t i = 0; i < ba_size; i++)
    {
        error_positions[i] = alice_bit_array[i] ^ bob_bit_array[i];
    }
    return error_positions;
}

void shuffle_bits(int* alice_bit_array, int* bob_bit_array, size_t ba_size, int seed) {
    mt19937 rng(seed);
    shuffle(alice_bit_array, alice_bit_array + ba_size, rng);
    rng.seed(seed);
    shuffle(bob_bit_array, bob_bit_array + ba_size, rng);
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

int* calculate_syndrome(int* bit_block, size_t block_size, umat& hash_matrix) {
    ucolvec column_vector(block_size);
    for (size_t i = 0; i < block_size; i++)
    {
        column_vector(i) = bit_block[i];
    }
    ucolvec syndrome = hash_matrix * column_vector;
 
    int* syndrome_array = new int [syndrome.size()];
    for (size_t i = 0; i < syndrome.size(); i++)
    {
        syndrome_array[i] = syndrome(i) % 2;
    }
    return syndrome_array;
}

void correct_error(int* bit_block, size_t block_size, int* first_syndrome, int* second_syndrome, size_t syndrome_power) {
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
    int* random_bit_array = generate_random_bit_array(SIFTED_KEY_LENGTH);
    _PRINT_ARRAY(random_bit_array, SIFTED_KEY_LENGTH, 8);
    int* bit_array_with_errors = introduce_errors(random_bit_array, SIFTED_KEY_LENGTH, 0.05);
    _PRINT_ARRAY(bit_array_with_errors, SIFTED_KEY_LENGTH, 8);
    //int* error_positions = calculate_error_positions(random_bit_array, bit_array_with_errors, SIFTED_KEY_LENGTH);
    
    /*int amp1_size = SIFTED_KEY_LENGTH - (int)(SIFTED_KEY_LENGTH / pow(2, INITIAL_SYNDROME_POWER));
    int* amp1_bit_array = new int[amp1_size];*/
    
    //auto syn_a = calculate_syndrome(random_bit_array, SIFTED_KEY_LENGTH, hash_matrix);
    //random_bit_array[0] = !random_bit_array[0];
    //print_array(random_bit_array, SIFTED_KEY_LENGTH, 7);
    /*auto syn_b = calculate_syndrome(random_bit_array, SIFTED_KEY_LENGTH, hash_matrix);
    print_array(syn_a, 3, 100);
    print_array(syn_b, 3, 100);
    correct_error(random_bit_array, SIFTED_KEY_LENGTH, syn_a, syn_b, 3);*/
    //print_array(random_bit_array, SIFTED_KEY_LENGTH, 7);
    size_t init_array_size = SIFTED_KEY_LENGTH;
    int syndrome_power = 3;
    size_t block_length = pow(2, syndrome_power);
    int blocks_count = SIFTED_KEY_LENGTH / block_length;
    umat hash_matrix = calculate_Hamming_hash_function(syndrome_power);

    vector<int> b_odd_parity_block_numbers;
    b_odd_parity_block_numbers.reserve(blocks_count);
    for (size_t i = 0; i < init_array_size; i += block_length)
    {
        if (block_parity_check(random_bit_array + i, block_length) != block_parity_check(bit_array_with_errors + i, block_length))
        {
            b_odd_parity_block_numbers.push_back(i / block_length);
        }
    }
    for (size_t i = 0; i < b_odd_parity_block_numbers.size(); i++)
    {
        cout << b_odd_parity_block_numbers[i]<<' ';
    }
    cout << endl;

    size_t amp1_array_size = init_array_size - (int)init_array_size / block_length;
    int* alice_privacy_amp_1 = new int[amp1_array_size];
    int* bob_privacy_amp_1 = new int[amp1_array_size];
    discard_bits_for_parity_check(random_bit_array, init_array_size, alice_privacy_amp_1, amp1_array_size, syndrome_power);
    discard_bits_for_parity_check(bit_array_with_errors, init_array_size, bob_privacy_amp_1, amp1_array_size, syndrome_power);
    _PRINT_ARRAY(alice_privacy_amp_1, amp1_array_size, 7);
    _PRINT_ARRAY(bob_privacy_amp_1, amp1_array_size, 7);

    block_length -= 1;
    int* alice_syndrome = new int[syndrome_power];
    int* bob_syndrome = new int[syndrome_power];
    for (size_t i = 0; i < b_odd_parity_block_numbers.size(); i++)
    {
        alice_syndrome = calculate_syndrome(alice_privacy_amp_1 + (b_odd_parity_block_numbers[i] * block_length), block_length, hash_matrix);
        bob_syndrome = calculate_syndrome(bob_privacy_amp_1 + (b_odd_parity_block_numbers[i] * block_length), block_length, hash_matrix);
        correct_error(alice_privacy_amp_1 + (b_odd_parity_block_numbers[i] * block_length), block_length, alice_syndrome, bob_syndrome, syndrome_power);
    }
    _PRINT_ARRAY(alice_privacy_amp_1, amp1_array_size, 7);
    _PRINT_ARRAY(bob_privacy_amp_1, amp1_array_size, 7);

    
    
    int remaining_bits_count = block_length - syndrome_power;
    size_t amp2_array_size = amp1_array_size - b_odd_parity_block_numbers.size() * syndrome_power;
    int* alice_privacy_amp_2 = new int[amp2_array_size];
    int* bob_privacy_amp_2 = new int[amp2_array_size];

    vector<int> discarded_bit_positions(syndrome_power + 1);
    for (size_t i = 0; i < discarded_bit_positions.size(); i++)
    {
        discarded_bit_positions[i] = pow(2, i) - 1;
    }
    int k = 0;
    int m = 0;
    for (size_t i = 0; i < b_odd_parity_block_numbers.size(); i++)
    {
        
    }

    //discard_bit_for_syndrome(bob_privacy_amp_1, amp1_array_size, test, discarded_bit_positions);
    
    /*
    auto start_time = steady_clock::now();
    for (size_t m = 0; m < 100000; m++)
    {
        discard_bits_for_parity_check(random_bit_array, SIFTED_KEY_LENGTH, amp1_bit_array, amp1_size, INITIAL_SYNDROME_POWER);
    }
    auto end_time = steady_clock::now();
    duration<double> execution_time = end_time - start_time;
    cout << "\n" << execution_time << endl;*/

   /* for (size_t i = 0; i < SIFTED_KEY_LENGTH; i++)
    {
        cout << random_bit_array[i];
    }
    cout << endl;
    for (size_t i = 0; i < SIFTED_KEY_LENGTH; i++)
    {
        cout << bit_array_with_errors[i];
    }
    cout << endl;
    for (size_t i = 0; i < SIFTED_KEY_LENGTH; i++)
    {
        cout << error_positions[i];
    }*/
    
    
    return 0;
}
