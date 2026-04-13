#include <omp.h>
#include <iostream>
#include <vector>
#include <memory>
#include <iomanip>
#include "state.hpp"
#include "model.hpp"
#include "lib/integrators.hpp"
#include "lib/solvers.hpp"
#include "lib/io.hpp"
#include "lib/modules.hpp"
#include "lib/config_reader.hpp"

using namespace num;
using namespace mod;
using namespace top;
using namespace mod::heat;

int main(int argc, char** argv) {
    std::string config_file = "input/heat_3d_implicit.txt";
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

    // 3D Grid dimensions
    int nx = config.get("nx", 21);
    int ny = config.get("ny", 21);
    int nz = config.get("nz", 21);
    double dx = config.get("dx", 0.5);
    double dy = config.get("dy", 0.5);
    double dz = config.get("dz", 0.5);
    double alpha = config.get("alpha", 0.1);
    double dt = config.get("dt", 0.1);
    double t_end = config.get("t_end", 2.0);

    Spatial3D spatial(nx, ny, nz, dx, dy, dz);
    auto state = std::make_shared<Heat3DImplicitState>(spatial, 20.0);
    
    // Set a hot "cube" in the center
    for (int k = 8; k <= 12; ++k) {
        for (int j = 8; j <= 12; ++j) {
            for (int i = 8; i <= 12; ++i) {
                state->temperatures[state->spatial.idx(i, j, k)] = 100.0;
            }
        }
    }

    // Boundary conditions (all faces at 20.0)
    auto model = std::make_shared<Heat3DModel>(alpha, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0);
    auto solver = std::make_shared<ConjugateGradientSolver>();
    auto integrator = std::make_shared<ImplicitEulerIntegrator>();

    StandardSimulator sim(model, state, solver, integrator);

    std::cout << "Starting 3D Implicit Heat Simulation\n";
    if (enable_vtk) std::cout << "VTK Export: Enabled\n";

    auto slicer_logger = [&spatial, enable_vtk](double t, const IState& s) {
        static int step_count = 0;
        const auto& h_state = dynamic_cast<const Heat3DImplicitState&>(s);
        
        if (enable_vtk) {
            char filename[100];
            std::sprintf(filename, "exports/heat_3d_%03d.vti", step_count);
            VTKExporter::export_vti_3d(filename, h_state.temperatures, (int)spatial.nx, (int)spatial.ny, (int)spatial.nz, 
                                       spatial.dx, spatial.dy, spatial.dz, "Temperature");
        }

        if (step_count % 5 == 0) {
            std::cout << "Time: " << std::fixed << std::setprecision(2) << t << " | ";
            
            // Log central slices cross-section temperature
            int mid_k = (int)h_state.spatial.nz / 2;
            int mid_j = (int)h_state.spatial.ny / 2;
            int mid_i = (int)h_state.spatial.nx / 2;
            double center_T = h_state.temperatures[h_state.spatial.idx(mid_i, mid_j, mid_k)];
            std::cout << "Center Temp: " << center_T << "\n";

            // Visual slice (XY slice at mid-Z)
            int skip = 2; // skip for visibility
            for (int j = 0; j < (int)h_state.spatial.ny; j += skip) {
                for (int i = 0; i < (int)h_state.spatial.nx; i += skip) {
                    double T = h_state.temperatures[h_state.spatial.idx(i, j, mid_k)];
                    if (T > 80.0) std::cout << "##";
                    else if (T > 50.0) std::cout << "++";
                    else if (T > 25.0) std::cout << "..";
                    else std::cout << "  ";
                }
                std::cout << "\n";
            }
            std::cout << "--------------------------------\n";
        }
        step_count++;
    };

    sim.run(t_end, dt, slicer_logger);

    std::cout << "3D Simulation Successful.\n";
    return 0;
}
