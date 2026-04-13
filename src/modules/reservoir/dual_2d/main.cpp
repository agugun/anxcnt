#include <iostream>
#include <vector>
#include <memory>
#include <iomanip>
#include "state.hpp"
#include "model.hpp"
#include "lib/integrators.hpp"
#include "lib/solvers.hpp"
#include "lib/io.hpp"
#include "modules/reservoir/well.hpp"
#include "modules/reservoir/reservoir_integrators.hpp"
#include "lib/config_reader.hpp"

using namespace num;
using namespace mod;
using namespace top;
using namespace mod::reservoir;

int main(int argc, char** argv) {
    std::string config_file = "input/reservoir_dual_2d.txt";
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] != '-') {
            config_file = argv[i];
            break;
        }
    }

    ConfigReader config;
    config.load(config_file);

    bool enable_vtk = config.get("enable_vtk", 0) != 0;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--vtk") enable_vtk = true;
    }

    // Reservoir Geometry
    int nx = config.get("nx", 20);
    int ny = config.get("ny", 20);
    double total_x = config.get("total_x", 2000.0);
    double total_y = config.get("total_y", 2000.0);
    double dx = total_x / nx, dy = total_y / ny;
    double h = config.get("h", 50.0);

    // Properties
    double k = config.get("k", 200.0);
    double phi = config.get("phi", 0.2);
    double mu_w = config.get("mu_w", 1.0);
    double mu_o = config.get("mu_o", 5.0);
    double initial_p = config.get("initial_p", 1000.0);
    double sw_res = config.get("sw_res", 0.2);

    Spatial2D spatial(nx, ny, dx, dy);
    auto state = std::make_shared<ReservoirDualPhase2DState>(spatial, initial_p, sw_res);

    // RelPerm helper for well instantiation
    auto rel_perm_func = [](double sw, double& krw, double& kro) {
        double sw_res = 0.2, so_res = 0.2;
        double swe = (sw - sw_res) / (1.0 - sw_res - so_res);
        swe = std::max(0.0, std::min(1.0, swe));
        krw = swe * swe;
        kro = (1.0 - swe) * (1.0 - swe);
    };

    auto idx_func = [nx](int i, int j) { return j * nx + i; };
    auto sw_func = [state, idx_func](int i, int j) { return state->water_saturations[idx_func(i, j)]; };

    // Wells: Waterflood pattern (ISourceSink Abstraction)
    std::vector<std::shared_ptr<ISourceSink>> sources;
    // Injector (Top-Left)
    sources.push_back(std::make_shared<ReservoirWellDual2D>(0, 0, 500.0, true, 
                                                           rel_perm_func, idx_func, sw_func, mu_w, mu_o));
    // Producer (Bottom-Right)
    sources.push_back(std::make_shared<ReservoirWellDual2D>(nx - 1, ny - 1, 500.0, false, 
                                                           rel_perm_func, idx_func, sw_func, mu_w, mu_o));

    auto model = std::make_shared<ReservoirDual2DModel>(k, phi, mu_w, mu_o, h, sources);

    auto solver = std::make_shared<ConjugateGradientSolver>();
    auto integrator = std::make_shared<ReservoirIMPESIntegrator<ReservoirDualPhase2DState, ReservoirDual2DModel>>();

    StandardSimulator sim(model, state, solver, integrator);

    std::cout << "Starting 2D Dual-Phase Reservoir Simulation\n";
    if (enable_vtk) std::cout << "VTK Export: Enabled\n";

    auto vtk_logger = [nx, ny, dx, dy, enable_vtk](double t, const IState& s) {
        static int step_count = 0;
        const auto& dual_state = dynamic_cast<const ReservoirDualPhase2DState&>(s);
        
        if (enable_vtk && step_count % 5 == 0) {
            char filename[100];
            std::sprintf(filename, "exports/reservoir_dual_%03d.vti", step_count / 5);
            
            std::vector<VTKField> fields = {
                {"Pressure", dual_state.pressures},
                {"WaterSaturation", dual_state.water_saturations}
            };
            VTKExporter::export_vti_multi_2d(filename, nx, ny, dx, dy, fields);
        }

        if (step_count % 5 == 0) {
            std::cout << "Time: " << std::fixed << std::setprecision(1) << t 
                      << " hr | Sw(Injector): " << dual_state.water_saturations[0] 
                      << " | Sw(Producer): " << dual_state.water_saturations[nx*ny-1] << "\n";
        }
        step_count++;
    };

    // Simulation for 200 days
    double dt = 12.0;      // 12 hours
    double t_end = 4800.0;  // 200 days
    sim.run(t_end, dt, vtk_logger);

    std::cout << "Simulation Successful.\n";
    return 0;
}
