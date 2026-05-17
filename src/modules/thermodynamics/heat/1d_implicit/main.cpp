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
    // 1. Controller: Configuration & Environment
    std::string config_file = "input/heat_1d_implicit.txt";
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] != '-') {
            config_file = argv[i];
            break;
        }
    }

    ConfigReader config;
    config.load(config_file);

    // Command line overrides
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--vtk") config.set("enable_vtk", 1);
        if (arg == "--csv") config.set("enable_csv", 1);
    }

    omp_set_num_threads(config.get("num_threads", 1));

    // 2. Initialize Specialized Simulation
    HeatSimulation sim;
    sim.build(config);

    // 3. Initial State
    auto st = sim.create_initial_state(config);

    // 4. Logger / Observer setup
    auto logger = std::make_shared<StandardLogger>(config);
    logger->set_grid(config.get("nx", 100), 1, 1, config.get("dx", 0.01));
    logger->add_field("Temperature", [](const IState& s) {
        return s.to_vector();
    });
    sim.add_observer(logger);

    // 5. Execution: Orchestration
    double t_end = config.get("t_end", 2.0);
    double dt = config.get("dt", 0.01);

    std::cout << "Starting Heat 1D Implicit Simulation\n";
    sim.run(t_end, dt, std::move(st));

    std::cout << "Simulation Successful.\n";
    return 0;
}
