#include <omp.h>
#include <iostream>
#include "simulation.hpp"
#include "lib/utils/config_reader.hpp"

using namespace mod::burgers;
using namespace utl;

int main(int argc, char** argv) {
    std::string config_file = "input/burgers_fdm.txt";
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] != '-') {
            config_file = argv[i];
            break;
        }
    }

    ConfigReader config;
    if (!config.load(config_file)) {
        std::cerr << "Failed to load config: " << config_file << "\n";
        return 1;
    }
    
    int num_threads = config.get("num_threads", 4);
    omp_set_num_threads(num_threads);

    auto [engine, state, logger] = BurgersSimulationBuilder::build(config);

    double dt = config.get("dt", 0.01);
    double t_end = config.get("t_end", 2.0);

    std::cout << "Starting 1D Burgers' Equation Simulation (FDM) [Refactored Architecture]...\n";
    engine->simulate(t_end, dt, std::move(state));

    std::cout << "Burgers Simulation Successfully Completed.\n";
    return 0;
}
