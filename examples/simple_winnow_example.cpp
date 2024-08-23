// #include <iostream>

// int main()
// {
//     int *a_array = new int[32]{0, 1, 0, 1, 1, 1, 0, 0,
//                                1, 1, 0, 1, 0, 1, 1, 1,
//                                0, 1, 0, 1, 1, 1, 0, 0,
//                                0, 1, 1, 1, 0, 1, 0, 1};
//     int *b_array = new int[32]{0, 1, 0, 1, 1, 1, 0, 0,
//                                1, 1, 0, 1, 0, 1, 1, 1,
//                                0, 1, 0, 1, 1, 1, 0, 0,
//                                0, 1, 1, 1, 0, 1, 0, 1};

//     b_array[10] = 1; // block #2 (one error)

//     b_array[17] = 0; // block #3 (two errors)
//     b_array[23] = 1;

//     b_array[24] = 1; // block #4 (three errors)
//     b_array[25] = 0;
//     b_array[26] = 0;

//     print_array(a_array, 32, 16);
//     print_array(b_array, 32, 16);

//     int *a_out_array = new int[32];
//     int *b_out_array = new int[32];
//     size_t out_len = winnow(a_array, b_array, 32, INITIAL_SYNDROME_POWER + 1, a_out_array, b_out_array);

//     print_array(a_out_array, out_len, 100);
//     print_array(b_out_array, out_len, 100);

//     delete[] a_array;
//     delete[] b_array;
//     delete[] a_out_array;
//     delete[] b_out_array;
// }