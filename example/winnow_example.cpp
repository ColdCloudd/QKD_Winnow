#include <iostream>

#include "config.hpp"
#include "utils.hpp"
#include "winnow_algorithm.hpp"

config_data CFG {};

int main()
{
    CFG.TRACE_WINNOW = true;

    const size_t ARR_LEN = 32;
    int *a_array = new int[ARR_LEN]{0, 1, 0, 1, 1, 1, 0, 0,
                                    1, 1, 0, 1, 0, 1, 1, 1,
                                    0, 1, 0, 1, 1, 1, 0, 0,
                                    0, 1, 1, 1, 0, 1, 0, 1};
    int *b_array = new int[ARR_LEN]{0, 1, 0, 1, 1, 1, 0, 0,
                                    1, 1, 0, 1, 0, 1, 1, 1,
                                    0, 1, 0, 1, 1, 1, 0, 0,
                                    0, 1, 1, 1, 0, 1, 0, 1};

    b_array[10] = b_array[10] ^ 1; // block #2 (one error)

    b_array[17] = b_array[17] ^ 1; // block #3 (two errors)
    b_array[23] = b_array[23] ^ 1;

    b_array[24] = b_array[24] ^ 1; // block #4 (three errors)
    b_array[25] = b_array[25] ^ 1;
    b_array[26] = b_array[26] ^ 1;

    int *a_out_array = new int[ARR_LEN];
    int *b_out_array = new int[ARR_LEN];
    size_t syndrome_power = 3;
    int **hash_mat = calculate_Hamming_hash_matrix(syndrome_power);

    winnow_result result = winnow(a_array, b_array, ARR_LEN, syndrome_power, hash_mat, a_out_array, b_out_array);

    for (size_t i = 0; i < syndrome_power; ++i)
    {
        delete[] hash_mat[i];
    }
    delete[] hash_mat;
    delete[] a_array;
    delete[] b_array;
    delete[] a_out_array;
    delete[] b_out_array;
}