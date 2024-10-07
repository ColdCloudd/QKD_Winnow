#include "bit_array_operations.hpp"

// Generates Alice's key
void generate_random_bit_array(std::mt19937 &prng, size_t length, int *const output_random_bit_array)
{
    std::uniform_int_distribution<int> distribution(0, 1);

    // Generate random bits and fill the vector
    for (int i = 0; i < length; ++i)
    {
        output_random_bit_array[i] = distribution(prng);
    }
}

// Generates Bob's key by making errors in Alice's key with a given QBER probability (Uniform distribution)
void introduce_errors(std::mt19937 &prng, const int *const bit_array, size_t array_length, float QBER,
                      int *const output_bit_array_with_errors)
{
    size_t num_errors = static_cast<size_t>(array_length * QBER);
    if (num_errors == 0)
    {
        std::copy(bit_array, bit_array + array_length, output_bit_array_with_errors);
    }
    else
    {
        size_t *error_positions = new size_t[array_length];
        for (size_t i = 0; i < array_length; ++i)
        {
            error_positions[i] = i;
        }

        std::shuffle(error_positions, error_positions + array_length, prng);
        std::copy(bit_array, bit_array + array_length, output_bit_array_with_errors);

        for (size_t i = 0; i < num_errors; ++i)
        {
            output_bit_array_with_errors[error_positions[i]] ^= 1;
        }

        delete[] error_positions;
    }
}

// Shuffles Alice's and Bob's bit arrays by seed
void shuffle_array_bits(int *const alice_bit_array, int *const bob_bit_array, size_t array_length, int seed)
{
    std::mt19937 rng(seed);
    std::shuffle(alice_bit_array, alice_bit_array + array_length, rng);
    rng.seed(seed);
    std::shuffle(bob_bit_array, bob_bit_array + array_length, rng);
}

// Calculates an array that consists of 0 and 1, where 1 denotes that Alice's and Bob's bits are different
void calculate_error_positions(const int *const alice_bit_array, const int *const bob_bit_array, size_t array_length,
                               int *const output_error_positions)
{
    for (size_t i = 0; i < array_length; i++)
    {
        output_error_positions[i] = alice_bit_array[i] ^ bob_bit_array[i];
    }
}