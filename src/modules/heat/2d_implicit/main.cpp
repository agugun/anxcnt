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
    std::string config_file = "input/heat_2d_implicit.txt";
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

    int nx = config.get("nx", 61);
    int ny = config.get("ny", 61);
    double dx = config.get("dx", 0.1);
    double dy = config.get("dy", 0.1);
    double alpha = config.get("alpha", 0.1);
    double dt = config.get("dt", 0.1);
    double t_end = config.get("t_end", 2.0);

    Spatial2D spatial(nx, ny, dx, dy);
    auto state = std::make_shared<Heat2DImplicitState>(spatial, 20.0);
    for (int j = 20; j <= 40; ++j) {
        for (int i = 20; i <= 40; ++i) {
            state->temperatures[state->spatial.idx(i, j)] = 100.0;
        }
    }

    auto model = std::make_shared<Heat2DModel>(alpha, 20.0, 20.0, 20.0, 20.0);
    auto solver = std::make_shared<ConjugateGradientSolver>();
    auto integrator = std::make_shared<ImplicitEulerIntegrator>();

    StandardSimulator sim(model, state, solver, integrator);

    std::cout << "Starting 2D Implicit Heat Simulation\n";
    if (enable_vtk) std::cout << "VTK Export: Enabled\n";

    auto grid_logger = [&spatial, enable_vtk](double t, const IState& s) {
        static int step_count = 0;
        const auto& h_state = dynamic_cast<const Heat2DImplicitState&>(s);

        if (enable_vtk) {
            char filename[100];
            std::sprintf(filename, "exports/heat_2d_%03d.vti", step_count);
            VTKExporter::export_vti_2d(filename, h_state.temperatures, (int)spatial.nx, (int)spatial.ny, spatial.dx, spatial.dy, "Temperature");
        }

        if (step_count % 5 == 0) {
            std::cout << "Time: " << std::fixed << std::setprecision(2) << t << "\n";
            
            // Print 10x10 subsample
            int rows = 10, cols = 10;
            int step_x = (int)h_state.spatial.nx / cols;
            int step_y = (int)h_state.spatial.ny / rows;

            for (int j = 0; j < (int)h_state.spatial.ny; j += step_y) {
                for (int i = 0; i < (int)h_state.spatial.nx; i += step_x) {
                    double T = h_state.temperatures[h_state.spatial.idx(i, j)];
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

    sim.run(t_end, dt, grid_logger);

    std::cout << "Simulation Successful.\n";
    return 0;
}
