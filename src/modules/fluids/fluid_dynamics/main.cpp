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
    std::string config_file = "input/fluid_fem.txt";
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
    FluidSimulation sim;
    sim.build(config);

    auto st = sim.create_initial_state(config);

    // 2. Logger / Observer setup
    auto logger = std::make_shared<StandardLogger>(config);
    logger->set_grid(config.get("nx", 10), config.get("ny", 10), 1, config.get("L", 1.0)/config.get("nx", 10));
    logger->add_field("U-Velocity", [](const IState& s) {
        const auto& fs = static_cast<const FluidState&>(s);
        return fs.u;
    });
    logger->add_field("P-Field", [](const IState& s) {
        const auto& fs = static_cast<const FluidState&>(s);
        return fs.p;
    });
    sim.add_observer(logger);

    // 3. Execution: Orchestration
    double dt = config.get("dt", 0.1);
    double t_end = config.get("t_end", 5.0);

    std::cout << "Starting 2D Fluid Dynamics Simulation (FEM)\n";
    sim.run(t_end, dt, std::move(st));

    std::cout << "Fluid Simulation Successfully Completed.\n";
    return 0;
}
