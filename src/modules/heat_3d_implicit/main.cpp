#include <iostream>
#include <vector>
#include <memory>
#include <iomanip>
#include "state.hpp"
#include "model.hpp"
#include "lib/solvers.hpp"
#include "lib/integrators.hpp"
#include "lib/io.hpp"

using namespace num;
using namespace mod;
using namespace top;
using namespace mod::physics_heat;

int main() {
    // 3D Grid dimensions
    int nx = 21, ny = 21, nz = 21;
    double dx = 0.5, dy = 0.5, dz = 0.5;
    double alpha = 0.1, dt = 0.1, t_end = 2.0;

    auto state = std::make_shared<Heat3DState>(nx, ny, nz, dx, dy, dz, 20.0);
    
    // Set a hot "cube" in the center
    for (int k = 8; k <= 12; ++k) {
        for (int j = 8; j <= 12; ++j) {
            for (int i = 8; i <= 12; ++i) {
                state->temperatures[state->idx(i, j, k)] = 100.0;
            }
        }
    }

    // Boundary conditions (all faces at 20.0)
    auto model = std::make_shared<Heat3DModel>(alpha, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0);
    auto solver = std::make_shared<ConjugateGradientSolver>();
    auto integrator = std::make_shared<ImplicitEulerIntegrator>();

    StandardSimulator sim(model, state, solver, integrator);

    std::cout << "Starting 3D Implicit Heat Simulation (Generic Simulator)\n";
    std::cout << "Grid: " << nx << "x" << ny << "x" << nz << " nodes\n";
    std::cout << "Exporting results to exports/heat_3d_*.vtk\n";

    auto slicer_logger = [nx, ny, nz, dx, dy, dz](double t, const IState& s) {
        static int step_count = 0;
        const auto& h_state = dynamic_cast<const Heat3DState&>(s);
        
        // Export VTI every step
        char filename[100];
        std::sprintf(filename, "exports/heat_3d_%03d.vti", step_count);
        VTKExporter::export_vti_3d(filename, h_state.temperatures, nx, ny, nz, dx, dy, dz, "Temperature");

        if (step_count % 5 == 0) {
            std::cout << "Time: " << std::fixed << std::setprecision(2) << t << " | ";
            
            // Log central slices cross-section temperature
            int mid_k = h_state.nz / 2;
            int mid_j = h_state.ny / 2;
            int mid_i = h_state.nx / 2;
            double center_T = h_state.temperatures[h_state.idx(mid_i, mid_j, mid_k)];
            std::cout << "Center Temp: " << center_T << "\n";

            // Visual slice (XY slice at mid-Z)
            int skip = 2; // skip for visibility
            for (int j = 0; j < h_state.ny; j += skip) {
                for (int i = 0; i < h_state.nx; i += skip) {
                    double T = h_state.temperatures[h_state.idx(i, j, mid_k)];
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
