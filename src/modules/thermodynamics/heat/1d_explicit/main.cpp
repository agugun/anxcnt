#include <omp.h>
#include <iostream>
#include <memory>
#include "simulation.hpp"
#include "lib/utils/config_reader.hpp"

using namespace utl;
using namespace mod::heat;

int main(int argc, char** argv) {
    // 1. Controller: Configuration & Environment
    std::string config_file = "input/heat_1d_explicit.txt";
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

    // 2. Construction: Simulation Engine & Logger
    auto [engine, state, logger] = Heat1DExplicitSimulation::build(config);

    // 3. Execution: Orchestration
    double t_end = config.get("t_end", 2.0);
    double dt = config.get("dt", 0.0001); // Small dt for explicit stability

    std::cout << "Starting Heat 1D Explicit Simulation [Modular Engine Architecture]\n";
    engine->simulate(t_end, dt, std::move(state));

    std::cout << "Simulation Successful.\n";
    return 0;
}