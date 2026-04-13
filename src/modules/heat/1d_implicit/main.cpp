#include <iostream>
#include <vector>
#include <memory>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "state.hpp"
#include "model.hpp"
#include "lib/solvers.hpp"
#include "lib/integrators.hpp"
#include "lib/io.hpp"
#include "lib/config_reader.hpp"

using namespace num;
using namespace mod;
using namespace top;
using namespace mod::heat;

int main(int argc, char** argv) {
    std::string config_file = "input/heat_1d_implicit.txt";
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] != '-') {
            config_file = argv[i];
            break;
        }
    }

    ConfigReader config;
    config.load(config_file);

    bool enable_vtk = config.get("enable_vtk", 0) != 0;
    bool enable_csv = config.get("enable_csv", 0) != 0;
    
    // Command line overrides
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--vtk") enable_vtk = true;
        else if (arg == "--csv") enable_csv = true;
    }

    size_t nx = config.get("nx", 100);
    double dx = 1.0 / (nx - 1);
    double alpha = config.get("alpha", 0.1);
    double dt = config.get("dt", 0.01);
    double t_end = config.get("t_end", 2.0);

    Spatial1D spatial(nx, dx);
    auto state = std::make_shared<Heat1DImplicitState>(spatial, 0.0);
    for (int i = 0; i < nx; ++i) {
        state->temperatures[i] = std::exp(-std::pow(i * dx - 0.5, 2) / 0.02);
    }

    auto model = std::make_shared<Heat1DModel>(alpha, 0.0, 0.0);
    auto solver = std::make_shared<LinearTridiagonalSolver>();
    auto integrator = std::make_shared<ImplicitEulerIntegrator>();

    StandardSimulator sim(model, state, solver, integrator);

    std::cout << "Starting Implicit Heat Simulation\n";
    if (enable_vtk) std::cout << "VTK Export: Enabled\n";
    if (enable_csv) std::cout << "CSV Export: Enabled\n";

    std::ofstream csv;
    if (enable_csv) {
        csv.open("exports/heat_results.csv");
        csv << "time,index,x,temp\n";
    }

    auto logger = [&csv, &spatial, enable_vtk, enable_csv](double t, const IState& s) {
        static int step = 0;
        const auto& h_state = dynamic_cast<const Heat1DImplicitState&>(s);
        
        // Export every 5 steps
        if (step % 5 == 0) {
            if (enable_vtk) {
                char filename[100];
                std::sprintf(filename, "exports/heat_1d_im_%03d.vti", step / 5);
                VTKExporter::export_vti_1d(filename, h_state.temperatures, (int)spatial.nx, spatial.dx, "Temperature");
            }

            if (enable_csv && csv.is_open()) {
                for (size_t i = 0; i < h_state.temperatures.size(); ++i) {
                    csv << t << "," << i << "," << i * spatial.dx << "," << h_state.temperatures[i] << "\n";
                }
            }
        }

        if (step % 20 == 0) {
            std::cout << "Time: " << std::fixed << std::setprecision(3) << t 
                      << " | Center Temp: " << h_state.temperatures[h_state.temperatures.size() / 2] << "\n";
        }
        
        step++;
    };

    sim.run(t_end, dt, logger);

    if (enable_csv && csv.is_open()) csv.close();
    std::cout << "Simulation Successful.\n";
    return 0;
}
