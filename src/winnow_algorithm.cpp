#include "winnow_algorithm.hpp"

// Discards bits at the first position in the block for privacy amplification
void discard_bits_for_parity_check(const int *const source_bit_array, const size_t &source_array_length,
                                   int *const destination_bit_array, const size_t &syndrome_power)
{
    size_t source_block_size = static_cast<size_t>(pow(2, syndrome_power));
    size_t destination_block_size = source_block_size - 1;
    for (size_t i = 0, j = 0; i < source_array_length; i += source_block_size, j += destination_block_size)
    {
        std::copy(source_bit_array + i + 1, source_bit_array + i + source_block_size, destination_bit_array + j);
    }
}

// Discards bits at positions 2^n, where n=(0,1, ... syndrome_power - 1) in the block for privacy amplification
void discard_bits_for_syndrome(const int *const source_bit_block, int *const destination_bit_block,
                               const std::vector<int> &discarded_bit_positions)
{
    int destination_current_start = 0;
    // Copying bits that are between the positions described in discarded_bit_positions
    for (size_t i = 1; i < discarded_bit_positions.size() - 1; i++)
    {
        destination_current_start += discarded_bit_positions[i] - discarded_bit_positions[i - 1] - 1;
        std::copy(source_bit_block + discarded_bit_positions[i] + 1, source_bit_block + discarded_bit_positions[i + 1], destination_bit_block + destination_current_start);
    }
}

// Calculates the parity bit for a block
bool calculate_block_parity(const int *const bit_block, const size_t &block_length)
{
    return std::accumulate(bit_block, bit_block + block_length, 0) % 2 == 0;
}

// Computes a matrix based on the Hamming hash function
int **calculate_Hamming_hash_matrix(size_t syndrome_power)
{
    size_t block_length = static_cast<size_t>(pow(2, syndrome_power)) - 1;

    int **hash_matrix = new int *[syndrome_power];
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
void calculate_syndrome(const int *const bit_block, const size_t &syndrome_power, const size_t &block_length,
                        const int *const *hash_matrix, int *const output_syndrome)
{
    int xor_sum = 0;
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
void correct_error(int *const bit_block, const int *const first_syndrome, const int *const second_syndrome,
                   const size_t &syndrome_power)
{
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

winnow_result winnow(int *const alice_bit_array, int *const bob_bit_array, size_t array_length, size_t syndrome_power,
              const int* const* hash_mat, int *const output_alice_bit_array, int *const output_bob_bit_array)
{
    
    size_t block_len = static_cast<size_t>(pow(2, syndrome_power));
    size_t blocks_cnt = array_length / block_len;

    std::vector<int> diff_par_blocks; // Contains the numbers of blocks whose parity bits did not match for Alice and Bob
    diff_par_blocks.reserve(blocks_cnt);
    for (size_t i = 0; i < array_length; i += block_len)
    {
        if (calculate_block_parity(alice_bit_array + i, block_len) != calculate_block_parity(bob_bit_array + i, block_len))
        {
            diff_par_blocks.push_back(static_cast<int>(i / block_len));
        }
    }

    if(CFG.TRACE_WINNOW)
    {
        fmt::print(fg(fmt::color::blue), "______________WINNOW_TRACE______________\n");
        fmt::print(fg(fmt::color::blue), "{}      - Numbers of blocks with detected errors\n", fmt::join(diff_par_blocks, " "));
    }
    
    size_t priv_amp_arr_len = array_length - static_cast<size_t>(array_length / block_len);
    int *alice_priv_amp = new int[priv_amp_arr_len];
    int *bob_priv_amp = new int[priv_amp_arr_len];
    // Privacy amplification by discarding the first bit in each block.
    discard_bits_for_parity_check(alice_bit_array, array_length, alice_priv_amp, syndrome_power);
    discard_bits_for_parity_check(bob_bit_array, array_length, bob_priv_amp, syndrome_power);

    if(CFG.TRACE_WINNOW)
    {
        print_array(alice_bit_array, array_length, block_len);
        fmt::print(fg(fmt::color::blue), "      - Alice's key\n");
        print_array(bob_bit_array, array_length, block_len);
        fmt::print(fg(fmt::color::blue), "      - Bob's key\n");

        print_array(alice_priv_amp, priv_amp_arr_len, block_len - 1);
        fmt::print(fg(fmt::color::blue), "      - Alice's key after first privacy maintenance\n");
        print_array(bob_priv_amp, priv_amp_arr_len, block_len - 1);
        fmt::print(fg(fmt::color::blue), "      - Bob's key after first privacy maintenance\n");
    }

    block_len -= 1;
    // Alice and Bob syndromes
    int *alice_syn = new int[syndrome_power];
    int *bob_syn = new int[syndrome_power];
    for (size_t i = 0; i < diff_par_blocks.size(); i++)
    {
        // Calculation of syndromes for blocks with non-matching parity bits, followed by error correction
        calculate_syndrome(alice_priv_amp + (diff_par_blocks[i] * block_len), syndrome_power, block_len, hash_mat, alice_syn);
        calculate_syndrome(bob_priv_amp + (diff_par_blocks[i] * block_len), syndrome_power, block_len, hash_mat, bob_syn);
        correct_error(alice_priv_amp + (diff_par_blocks[i] * block_len), alice_syn, bob_syn, syndrome_power);
    }

    delete[] alice_syn;
    delete[] bob_syn;

    if(CFG.TRACE_WINNOW)
    {
        print_array(alice_priv_amp, priv_amp_arr_len, block_len);
        fmt::print(fg(fmt::color::blue), "      - Alice's key after error correction\n");
        print_array(bob_priv_amp, priv_amp_arr_len, block_len);
        fmt::print(fg(fmt::color::blue), "      - Bob's key after error correction\n");
    }

    size_t remain_bits_cnt = block_len - syndrome_power; // Number of remaining bits in blocks for which syndromes were calculated
    size_t out_arr_len = priv_amp_arr_len - diff_par_blocks.size() * syndrome_power;

    // Contains bounds that specify valid bits [from 2^0, 2^1, ... , 2^(syndrome_power-1)],
    // and the last element (2^syndrome_power) which is the right boundary
    std::vector<int> disc_bit_pos(syndrome_power + 1);
    for (size_t i = 0; i < disc_bit_pos.size(); i++)
    {
        disc_bit_pos[i] = (int)(pow(2, i) - 1);
    }

    size_t j = 0;            // Used to move through diff_par_blocks
    size_t copy_delta = 0;   // Specifies the number of bits that can be copied
    size_t dest_cur_pos = 0; // The current position in the destination array from which new bits can be inserted
    size_t diff_par_len = diff_par_blocks.size();
    for (size_t i = 0; i < priv_amp_arr_len;)
    {
        if (j < diff_par_len && (i / block_len == diff_par_blocks[j])) // In a block whose number is in diff_par_blocks, bits are discarded
        {
            discard_bits_for_syndrome(alice_priv_amp + i, output_alice_bit_array + dest_cur_pos, disc_bit_pos);
            discard_bits_for_syndrome(bob_priv_amp + i, output_bob_bit_array + dest_cur_pos, disc_bit_pos);
            dest_cur_pos += remain_bits_cnt;
            i += block_len;
            j++;
        }
        else if (diff_par_len == 0) // No errors in arrays -> copy the array completely
        {
            std::copy(alice_priv_amp, alice_priv_amp + priv_amp_arr_len, output_alice_bit_array);
            std::copy(bob_priv_amp, bob_priv_amp + priv_amp_arr_len, output_bob_bit_array);
            i += priv_amp_arr_len;
        }
        else if (j >= diff_par_len) // Remaining blocks with no errors are copied
        {
            copy_delta = priv_amp_arr_len - i;
            std::copy(alice_priv_amp + i, alice_priv_amp + i + copy_delta, output_alice_bit_array + dest_cur_pos);
            std::copy(bob_priv_amp + i, bob_priv_amp + i + copy_delta, output_bob_bit_array + dest_cur_pos);
            i += copy_delta;
        }
        else
        {
            copy_delta = diff_par_blocks[j] * block_len - i; // Blocks between two erroneous blocks are copied
            std::copy(alice_priv_amp + i, alice_priv_amp + i + copy_delta, output_alice_bit_array + dest_cur_pos);
            std::copy(bob_priv_amp + i, bob_priv_amp + i + copy_delta, output_bob_bit_array + dest_cur_pos);
            dest_cur_pos += copy_delta;
            i += copy_delta;
        }
    }

    // Used to account for the number of bits to be removed for padding
    bool last_block_with_error = false; 
    if(diff_par_len > 0 && (diff_par_blocks.back() == (blocks_cnt - 1)))
    {
        last_block_with_error = true;
    }

    if(CFG.TRACE_WINNOW)
    {
        print_array(output_alice_bit_array, out_arr_len, block_len);
        fmt::print(fg(fmt::color::blue), "      - Alice's key after second privacy maintenance\n");
        print_array(output_bob_bit_array, out_arr_len, block_len);
        fmt::print(fg(fmt::color::blue), "      - Bob's key after second privacy maintenance\n");
        fmt::print(fg(fmt::color::blue), "______________WINNOW_TRACE______________\n");
    }

    delete[] alice_priv_amp;
    delete[] bob_priv_amp;

    return {out_arr_len, last_block_with_error};
}