#include <omp.h>
#include <iostream>
#include <vector>
#include <memory>
#include <iomanip>
#include "state.hpp"
#include "model.hpp"
#include "lib/solvers.hpp"
#include "lib/integrators.hpp"
#include "lib/io.hpp"
#include "modules/reservoir/well.hpp"
#include "lib/config_reader.hpp"

using namespace num;
using namespace mod;
using namespace top;
using namespace mod::reservoir;

int main(int argc, char** argv) {
    std::string config_file = "input/reservoir_3d.txt";
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

    bool enable_vtk = config.get("enable_vtk", 0) != 0;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--vtk") enable_vtk = true;
    }

    // Reservoir Parameters
    int nx = config.get("nx", 21);
    int ny = config.get("ny", 21);
    int nz = config.get("nz", 11);
    double total_x = config.get("total_x", 2100.0);
    double total_y = config.get("total_y", 2100.0);
    double total_z = config.get("total_z", 110.0);
    double dx = total_x / nx, dy = total_y / ny, dz = total_z / nz;
    
    double k = config.get("k", 100.0);
    double phi = config.get("phi", 0.2);
    double mu = config.get("mu", 1.0);
    double ct = config.get("ct", 1e-5);
    double B = config.get("B", 1.0);
    double initial_p = config.get("initial_p", 5000.0);
    double q_well = config.get("q_well", 500.0);
    double dt = config.get("dt", 24.0);
    double t_end = config.get("t_end", 720.0);

    // Well term conversion: (5.615 / 24.0) * B / (phi * ct * dx * dy * dz)
    double q_scale = (5.615 / 24.0) * B / (phi * dx * dy * dz * ct);
    auto idx_func = [nx, ny](int i, int j, int k) { return k * (nx * ny) + j * nx + i; };

    // Central producer completed across ALL layers (ISourceSink Abstraction)
    std::vector<std::shared_ptr<ISourceSink>> sources;
    sources.push_back(std::make_shared<ConstantRateWell>(nx / 2, ny / 2, 0, nz - 1, 500.0, q_scale, idx_func));

    // Create well location mask
    std::vector<double> well_mask(nx * ny * nz, 0.0);
    for (auto& s : sources) {
        auto w_ptr = std::dynamic_pointer_cast<IWell>(s);
        if (w_ptr) {
            for (int k = w_ptr->k_min; k <= w_ptr->k_max; ++k) {
                well_mask[k * (nx * ny) + w_ptr->j * nx + w_ptr->i] = 1.0;
            }
        }
    }

    // Simulation settings

    Spatial3D spatial(nx, ny, nz, dx, dy, dz);
    auto state = std::make_shared<Reservoir3DState>(spatial, initial_p);
    auto model = std::make_shared<Reservoir3DModel>(k, phi, mu, ct, B, sources);
    
    auto solver = std::make_shared<ConjugateGradientSolver>();
    auto integrator = std::make_shared<ImplicitEulerIntegrator>();

    StandardSimulator sim(model, state, solver, integrator);

    auto first_well = std::dynamic_pointer_cast<mod::IWell>(sources[0]);
    std::cout << "Starting 3D Reservoir Simulation\n";
    if (enable_vtk) std::cout << "VTK Export: Enabled\n";

    auto vtk_logger = [nx, ny, nz, dx, dy, dz, &well_mask, first_well, enable_vtk](double t, const IState& s) {
        static int step_count = 0;
        const auto& r_state = dynamic_cast<const Reservoir3DState&>(s);
        
        if (enable_vtk && step_count % 5 == 0) {
            char filename[100];
            std::sprintf(filename, "exports/reservoir_3d_%03d.vti", step_count / 5);
            
            std::vector<VTKField> fields = {
                {"Pressure", r_state.pressures},
                {"WellLocation", well_mask}
            };
            VTKExporter::export_vti_multi_3d(filename, nx, ny, nz, dx, dy, dz, fields);
        }

        if (step_count % 5 == 0) {
            int c_well = first_well->k_min * (nx * ny) + first_well->j * nx + first_well->i;
            std::cout << "Time: " << std::fixed << std::setprecision(1) << t 
                      << " hr | Well Pressure: " << r_state.pressures[c_well] << " psi\n";
        }
        step_count++;
    };

    sim.run(t_end, dt, vtk_logger);

    std::cout << "Simulation Successful.\n";
    return 0;
}
