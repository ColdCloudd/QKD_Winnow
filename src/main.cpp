#include "config.hpp"
#include "utils.hpp"
#include "simulation.hpp"

config_data CFG {};

int main()
{
    try
    {
        fs::path config_path =  fs::path(SOURCE_DIR) / "config.json";
        CFG = get_config_data(config_path);
        
        std::vector<std::vector<size_t>> trial_combinations {};
        if (CFG.USE_SPECIFIED_COMBINATIONS)
        {
            trial_combinations = CFG.COMBINATIONS;
        }
        else
        {
            trial_combinations = cartesian_product(CFG.COMBINATION_ELEMENTS);
        }
        std::vector<test_combination> combinations = prepare_combinations(trial_combinations, CFG.QBER);
        std::vector<test_result> result = run_simulation(combinations);

        fs::path result_dir_path = fs::path(SOURCE_DIR) / "results";
        if (!fs::exists(result_dir_path))
        {
            fs::create_directories(result_dir_path);
        }
        fmt::print(fg(fmt::color::green),"All tests were completed successfully! \nThe results will be written to the directory: {}\n", result_dir_path.string());
        write_file(result, result_dir_path);
    }
    catch(const std::exception& e)
    {
        fmt::print(stderr, fg(fmt::color::red), "ERROR: {}\n", e.what());
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
