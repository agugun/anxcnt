#include <iostream>
#include <vector>
#include <memory>
#include <omp.h>
#include "state.hpp"
#include "model.hpp"
#include "lib/solvers.hpp"
#include "lib/integrators.hpp"
#include "lib/io.hpp"
#include "lib/operators.hpp"
#include "lib/config_reader.hpp"

using namespace top;
using namespace mod::reservoir;
using namespace num;

int main(int argc, char** argv) {
    std::string config_file = "input/reservoir_oil_gas_2d.txt";
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] != '-') {
            config_file = argv[i];
            break;
        }
    }

    ConfigReader config;
    config.load(config_file);

    int num_threads = config.get("num_threads", 1);
    omp_set_num_threads(num_threads);
    std::cout << "OpenMP Active Threads: " << num_threads << " (Max Available: " << omp_get_max_threads() << ")" << std::endl;

    bool enable_vtk = config.get("enable_vtk", 0) != 0;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--vtk") enable_vtk = true;
    }

    std::cout << "--- 2D Dual-Phase Oil and Gas Simulation ---" << std::endl;
    if (enable_vtk) std::cout << "VTK Export: Enabled" << std::endl;

    // 1. Setup Grid
    int nx = config.get("nx", 20);
    int ny = config.get("ny", 20);
    double dx = config.get("dx", 100.0);
    double dy = config.get("dy", 100.0);
    Spatial2D spatial(nx, ny, dx, dy);
    double initial_sg = config.get("initial_sg", 0.15);
    auto state = std::make_shared<ReservoirOilGas2DState>(spatial, config.get("initial_p", 5000.0), initial_sg);

    // 2. Define Wells
    std::vector<std::shared_ptr<mod::ISourceSink>> wells;
    
    // Producer
    auto prod = std::make_shared<ReservoirWellOilGas2D>(
        nx-1, ny-1, config.get("q_rate", 200.0), false,
        [&](double sg, double& krog, double& krg) {
            // Simplified Corey for Gas-Oil
            double sge = sg / (1.0 - 0.2 - 0.1); 
            sge = std::max(0.0, std::min(1.0, sge));
            krg = sge * sge;
            krog = (1.0 - sge) * (1.0 - sge);
        },
        config.get("mu_o", 2.0), config.get("phi", 0.2), config.get("mu_g", 0.02),
        [&](double p){ return 1.0 * (14.7 / std::max(1.0, p)); }
    );
    
    wells.push_back(prod);

    // 3. Setup Model
    auto model = std::make_shared<ReservoirOilGas2DModel>(config.get("k", 100.0), config.get("phi", 0.2), config.get("mu_o", 2.0), config.get("mu_g", 0.02), config.get("h", 50.0), wells);

    // 4. Setup Simulator
    auto solver = std::make_shared<NewtonSolver>(); 
    auto integrator = std::make_shared<FullyImplicitIntegrator>();
    StandardSimulator sim(model, state, solver, integrator);

    // 5. Run Simulation
    double t_end = config.get("t_end", 2400.0);
    double dt = config.get("dt", 1.0); 

    sim.run(t_end, dt, [&](double t, const IState& s) {
        const auto& cur_state = dynamic_cast<const ReservoirOilGas2DState&>(s);
        
        // Console Telemetry
        if ((int)t % 24 == 0) {
            double avg_p = 0.0;
            double avg_sg = 0.0;
            for(double p : cur_state.pressures) avg_p += p;
            for(double sg : cur_state.gas_saturations) avg_sg += sg;
            avg_p /= cur_state.pressures.size();
            avg_sg /= cur_state.gas_saturations.size();

            std::cout << "Time: " << t << " hr | Avg P: " << avg_p << " psi | Avg Sg: " << avg_sg << std::endl;
        }

        // VTK Export
        if (enable_vtk && (int)t % 24 == 0) {
            std::string filename = "exports/reservoir_oil_gas_" + std::to_string((int)t) + ".vti";
            std::vector<VTKField> fields = {
                {"Pressure", cur_state.pressures},
                {"GasSaturation", cur_state.gas_saturations}
            };
            VTKExporter::export_vti_multi_2d(filename, cur_state.spatial.nx, cur_state.spatial.ny, cur_state.spatial.dx, cur_state.spatial.dy, fields);
        }
    });

    std::cout << "Simulation Successful." << std::endl;
    return 0;
}
