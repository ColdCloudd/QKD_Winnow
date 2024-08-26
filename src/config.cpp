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
            throw std::runtime_error("Number of threads must be greater than or equal to one!");
        }

        cfg.TRIALS_NUMBER = config["trials_number"].template get<size_t>();
        if (cfg.TRIALS_NUMBER < 1)
        {
            throw std::runtime_error("Number of trials must be greater than or equal to one!");
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

        cfg.INITIAL_SYNDROME_POWER = config["initial_syndrome_power"].template get<size_t>();
        if (cfg.INITIAL_SYNDROME_POWER < 3)
        {
            throw std::runtime_error("Minimum initial syndrome power is 3!");
        }
        
        cfg.QBER = config["qber"].template get<std::vector<double>>();
        for (size_t i = 0; i < cfg.QBER.size(); i++)
        {
            if (cfg.QBER[i] < 0 || cfg.QBER[i] > 1)
            {
                throw std::runtime_error("QBER cannot be less than zero or greater than one!");
            }
        }

        cfg.COMBINATION_ELEMENTS = config["combination_elements"].template get<std::vector<std::vector<size_t>>>();
        cfg.TRACE_WINNOW = config["trace_winnow"].template get<bool>();
        return cfg;
    }
    catch(const std::exception& e)
    {
        fmt::print(stderr, fg(fmt::color::red),"An error occurred while reading a configuration parameter.\n");
        throw;
    }
}