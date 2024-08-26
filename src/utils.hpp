#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <filesystem>

#include <fmt/core.h>
#include <fmt/color.h>
#include <fmt/ranges.h>

namespace fs = std::filesystem;

void print_array(const int *const bit_array, size_t array_length, size_t block_length);

