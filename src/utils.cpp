#include "utils.hpp"

void print_array(const int *const bit_array, size_t array_length, size_t block_length)
{
    for (size_t i = 0; i < array_length; i++)
    {
        if (i % block_length == 0 && i != 0)
        {
            fmt::print(" ");
        }
        fmt::print(fg(fmt::color::blue), "{}", bit_array[i]);
    }
}

