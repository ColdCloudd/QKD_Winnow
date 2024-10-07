#include "config.hpp"

config_data get_config_data(fs::path config_path)
{
    if (!fs::exists(config_path))
    {
        throw std::runtime_error("Configuration file not found: " + config_path.string());
    }

    std::ifstream config_file(config_path);
    if (!config_file.is_open())
    {
        throw std::runtime_error("Failed to open configuration file: " + config_path.string());
    }

    json config = json::parse(config_file);
    config_file.close();
    if (config.empty())
    {
        throw std::runtime_error("Configuration file is empty: " + config_path.string());
    }
    
    try
    {
        config_data cfg{};
        cfg.THREADS_NUMBER = config["threads_number"].template get<size_t>();
        if (cfg.THREADS_NUMBER < 1)
        {
            throw std::runtime_error("Number of threads must be >= 1!!");
        }

        cfg.TRIALS_NUMBER = config["trials_number"].template get<size_t>();
        if (cfg.TRIALS_NUMBER < 1)
        {
            throw std::runtime_error("Number of trials must be >= 1!!");
        }

        if (config["use_config_simulation_seed"].template get<bool>())
        {
            cfg.SIMULATION_SEED = config["simulation_seed"].template get<size_t>();
        }
        else
        {
            cfg.SIMULATION_SEED = time(nullptr);
        }
        
        cfg.SHUFFLE_MODE = config["shuffle_mode"].template get<bool>();

        cfg.SIFTED_KEY_LENGTH = config["sifted_key_length"].template get<size_t>();
        if (cfg.SIFTED_KEY_LENGTH < 8)
        {
            throw std::runtime_error("Minimum sifted key length is 8!");
        }

        cfg.INITIAL_SYNDROME_LENGTH = config["initial_syndrome_length"].template get<size_t>();
        if (cfg.INITIAL_SYNDROME_LENGTH < 3)
        {
            throw std::runtime_error("Minimum initial syndrome length is 3!");
        }
        
        cfg.QBER = config["qber"].template get<std::vector<double>>();
        for (size_t i = 0; i < cfg.QBER.size(); i++)
        {
            if (cfg.QBER[i] <= 0. || cfg.QBER[i] >= 1.)
            {
                throw std::runtime_error("QBER must be: 0 < QBER < 1!");
            }
        }

        cfg.USE_SPECIFIED_COMBINATIONS = config["use_specified_combinations"].template get<bool>();
        if (cfg.USE_SPECIFIED_COMBINATIONS)
        {
            cfg.COMBINATIONS = config["combinations"].template get<std::vector<std::vector<size_t>>>();
            if (cfg.COMBINATIONS.empty())
            {
                throw std::runtime_error("Vector 'combinations' cannot be empty!");
            }
            for (size_t i = 0; i < cfg.COMBINATIONS.size(); i++)
            {
                if (cfg.COMBINATIONS[i].empty())
                {
                    throw std::runtime_error("Vector 'combinations' cannot be empty!");
                }
            }
        }
        else
        {
            cfg.COMBINATION_ELEMENTS = config["combination_elements"].template get<std::vector<std::vector<size_t>>>();
            if (cfg.COMBINATION_ELEMENTS.empty())
            {
                throw std::runtime_error("Vector 'combination_elements' cannot be empty!");
            }
            for (size_t i = 0; i < cfg.COMBINATION_ELEMENTS.size(); i++)
            {
                if (cfg.COMBINATION_ELEMENTS[i].empty())
                {
                    throw std::runtime_error("Vector 'combination_elements' cannot be empty!");
                }
            }
        }
        
        cfg.TRACE_WINNOW = config["trace_winnow"].template get<bool>();
        return cfg;
    }
    catch(const std::exception& e)
    {
        fmt::print(stderr, fg(fmt::color::red),"An error occurred while reading a configuration parameter.\n");
        throw;
    }
}