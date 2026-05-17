#include <omp.h>
#include <iostream>
#include "simulation.hpp"
#include "lib/utils/config_reader.hpp"
#include "lib/utils/logger.hpp"

using namespace mod;
using namespace sim;
using namespace utl;
using namespace num;






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

    // 1. Construction: Simulation Engine & Logger
    BurgersSimulation sim;
    sim.build(config);

    auto st = sim.create_initial_state(config);

    // 2. Logger / Observer setup
    auto logger = std::make_shared<StandardLogger>(config);
    logger->set_grid(config.get("nx", 100), 1, 1, config.get("dx", 0.01));
    logger->add_field("Velocity", [](const IState& s) {
        return s.to_vector();
    });
    sim.add_observer(logger);

    // 3. Execution: Orchestration
    double dt = config.get("dt", 0.01);
    double t_end = config.get("t_end", 2.0);

    std::cout << "Starting 1D Burgers' Equation Simulation (FDM)\n";
    sim.run(t_end, dt, std::move(st));

    std::cout << "Burgers Simulation Successfully Completed.\n";
    return 0;
}
