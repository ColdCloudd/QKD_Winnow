#include "config.hpp"
#include "utils.hpp"
#include "simulation.hpp"

config_data CFG {};

const fs::path CONFIG_PATH = fs::path(SOURCE_DIR) / "config.json";
const fs::path RESULTS_DIRECTORY_PATH = fs::path(SOURCE_DIR) / "results";

int main()
{
    try
    {
        
        CFG = get_config_data(CONFIG_PATH);
        
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

        fmt::print(fg(fmt::color::green),"All tests were completed successfully! \nThe results will be written to the directory: {}\n", RESULTS_DIRECTORY_PATH.string());
        write_file(result, RESULTS_DIRECTORY_PATH);
    }
    catch(const std::exception& e)
    {
        fmt::print(stderr, fg(fmt::color::red), "ERROR: {}\n", e.what());
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
