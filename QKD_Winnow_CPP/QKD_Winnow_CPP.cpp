#include <cstdlib>
#include <iostream>
#include <numeric>
#include <vector>
#include <string>
#include <ctime>
#include <random>
#include <chrono>

using namespace std;
using namespace chrono;

int INITIAL_SYNDROME_POWER = 3;
int SIFTED_KEY_LENGTH = 50;
vector<float> BER = { 0.01, 0.03, 0.05, 0.06, 0.08, 0.1, 0.12, 0.15 };
vector<vector<int>> COMBINATION_ELEMENTS {{ 0, 1, 2 },
                                                { 0, 1, 2 },
                                                { 0, 1, 2 }};
int TRIAL_NUMBER = 10;
int PROCESSES_NUMBER = 12;
bool SHUFFLE_MODE = true;
string RESULT_FILE_PATH = "C:\\Users\\AGENT\\Desktop\\QKD_Winnow_NumPy\\results\\winnow_results.csv";


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

int* generate_random_bit_array(int length) {
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

int* introduce_errors(int* bit_array, int ba_size, float initial_error_probability) {
    mt19937 generator(time(0));
    uniform_real_distribution<double> distribution(0.0, 1.0);

    int* bit_array_with_errors = new int[ba_size];
    for (int i = 0; i < ba_size; ++i) {

        bit_array_with_errors[i] = bit_array[i] ^ (distribution(generator) < initial_error_probability);
    }

    return bit_array_with_errors;
}

void discard_bits_for_parity_check(int *source_bit_array, int sba_size, int *destination_bit_array, int dba_size, int syndrome_power) {
    int in_block_size = pow(2, syndrome_power);
    int out_block_size = in_block_size - 1;
    int j = 0;
    for (size_t i = 0, j = 0; i < sba_size && j < dba_size; i += in_block_size, j += out_block_size) {
        copy(source_bit_array + i + 1, source_bit_array + i + in_block_size, destination_bit_array + j);
    }
}

bool block_parity_check(int* bit_block, int block_size) {
    int sum = accumulate(bit_block, bit_block + block_size, 0);
    return sum % 2 == 0;
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
    /*block_parity_check(random_bit_array, 8);*/
    int* bit_array_with_errors = introduce_errors(random_bit_array, SIFTED_KEY_LENGTH, 0.05);
    
    /*int amp1_size = SIFTED_KEY_LENGTH - (int)(SIFTED_KEY_LENGTH / pow(2, INITIAL_SYNDROME_POWER));
    int* amp1_bit_array = new int[amp1_size];

    auto start_time = steady_clock::now();

    for (size_t m = 0; m < 100000; m++)
    {
        discard_bits_for_parity_check(random_bit_array, SIFTED_KEY_LENGTH, amp1_bit_array, amp1_size, INITIAL_SYNDROME_POWER);
    }
    auto end_time = steady_clock::now();
    duration<double> execution_time = end_time - start_time;
    cout << "\n" << execution_time << endl;*/

    for (size_t i = 0; i < SIFTED_KEY_LENGTH; i++)
    {
        cout << random_bit_array[i];
    }
    cout << endl;
    for (size_t i = 0; i < SIFTED_KEY_LENGTH; i++)
    {
        cout << bit_array_with_errors[i];
    }
    
    
    return 0;
}
