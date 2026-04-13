#include <iostream>
#include <vector>
#include <memory>
#include <fstream>
#include <iomanip>
#include "state.hpp"
#include "model.hpp"
#include "lib/solvers.hpp"
#include "lib/integrators.hpp"
#include "modules/reservoir/reservoir_integrators.hpp"
#include "lib/operators.hpp"
#include "modules/reservoir/well.hpp"
#include "lib/config_reader.hpp"

using namespace num;
using namespace mod;
using namespace top;
using namespace mod::reservoir;

int main(int argc, char** argv) {
    std::string config_file = "input/reservoir_black_oil_3d.txt";
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] != '-') {
            config_file = argv[i];
            break;
        }
    }

    ConfigReader config;
    config.load(config_file);

    std::cout << "--- Initiating 3D Black Oil Simulation ---\n";

    // 1. Grid Definition
    int nx = config.get("nx", 10);
    int ny = config.get("ny", 10);
    int nz = config.get("nz", 3);
    double dx = config.get("dx", 100.0);
    double dy = config.get("dy", 100.0);
    double dz = config.get("dz", 20.0);
    Spatial3D spatial(nx, ny, nz, dx, dy, dz);

    // 2. Initial State
    auto state = std::make_shared<ReservoirBlackOil3DState>(spatial);

    // 3. Wells Configuration
    auto rel_perm_wrapper = [](double sw, double sg, double& krw, double& kro, double& krg) {
        // Dummy relperm for well calc, should match model
        double swc = 0.2, sorg = 0.1, sgr = 0.05;
        double swe = (sw - swc) / (1.0 - swc - sorg);
        krw = std::pow(std::max(0.0, std::min(1.0, swe)), 2.0);
        double sge = (sg - sgr) / (1.0 - swc - sgr);
        krg = std::pow(std::max(0.0, std::min(1.0, sge)), 2.0);
        double so = 1.0 - sw - sg;
        double soe = (so - sorg) / (1.0 - swc - sorg);
        kro = std::pow(std::max(0.0, std::min(1.0, soe)), 3.0);
    };

    auto idx_wrapper = [&spatial](int i, int j, int k) { return spatial.idx(i, j, k); };
    auto var_accessor = [state](int i, int j, int k, int var) {
        int c = state->idx(i, j, k);
        if (var == 0) return state->p(c);
        if (var == 1) return state->sw(c);
        return state->sg(c);
    };

    std::vector<std::shared_ptr<ISourceSink>> sources;
    // Water Injector at (0,0) through all layers
    sources.push_back(std::make_shared<ReservoirWellBlackOil3D>(0, 0, 0, nz - 1, 500.0, true, 
                                                               rel_perm_wrapper, idx_wrapper, var_accessor));
    // Producer at (nx-1, ny-1) through all layers
    sources.push_back(std::make_shared<ReservoirWellBlackOil3D>(nx - 1, ny - 1, 0, nz - 1, 500.0, false, 
                                                               rel_perm_wrapper, idx_wrapper, var_accessor));

    // 4. Model & Solvers
    auto model = std::make_shared<ReservoirBlackOil3DModel>(100.0, 0.2, sources);
    auto solver = std::make_shared<NewtonSolver>();
    auto integrator = std::make_shared<FullyImplicitIntegrator>();

    StandardSimulator sim(model, state, solver, integrator);

    // 5. Simulation Run
    double dt = config.get("dt", 0.1);
    int steps = config.get("steps", 50);
    double t_end = steps * dt;

    std::cout << "Step | Avg Pressure | Avg Sw | Avg Sg | Time\n";
    std::cout << "-------------------------------------------\n";

    sim.run(t_end, dt, [&](double t, const IState& s) {
        static int current_step = 0;
        if (current_step == 0) { current_step++; return; } // Skip initial
        
        const auto& cur_state = dynamic_cast<const ReservoirBlackOil3DState&>(s);
        
        // Stat collection
        double sum_p = 0, sum_sw = 0, sum_sg = 0;
        int n_cells = nx * ny * nz;
        for (int i = 0; i < n_cells; ++i) {
            sum_p += cur_state.p(i);
            sum_sw += cur_state.sw(i);
            sum_sg += cur_state.sg(i);
        }

        std::cout << std::setfill(' ') << std::setw(4) << current_step << " | " 
                  << std::setw(12) << sum_p / n_cells << " | "
                  << std::setw(6) << sum_sw / n_cells << " | "
                  << std::setw(6) << sum_sg / n_cells << " | "
                  << "OK\n";
        current_step++;
    });

    std::cout << "\n3D Black Oil Simulation Completed Successfully.\n";
    return 0;
}
