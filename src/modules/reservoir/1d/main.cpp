#include <iostream>
#include <vector>
#include <memory>
#include <fstream>
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
    std::string config_file = "input/reservoir_1d.txt";
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
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--vtk") enable_vtk = true;
        else if (arg == "--csv") enable_csv = true;
    }

    // Reservoir Geometry
    int nx = config.get("nx", 50);
    double total_length = config.get("total_length", 2000.0);
    double dx = total_length / nx;
    double area = config.get("area", 10000.0);

    // Rock/Fluid properties
    double k = config.get("k", 100.0);
    double phi = config.get("phi", 0.2);
    double mu = config.get("mu", 1.0);
    double ct = config.get("ct", 1e-5);
    double B = config.get("B", 1.0);
    double initial_p = config.get("initial_p", 5000.0);

    // Well parameters
    int well_idx = nx / 2;
    double q_well = config.get("q_well", 200.0);
    // Simulation time
    double dt = config.get("dt", 12.0);
    double t_end = config.get("t_end", 720.0);

    // Well term conversion: (5.615 / 24.0) * B / (phi * ct * area * dx)
    double q_scale = (5.615 / 24.0) * B / (phi * ct * area * dx);
    auto idx_func = [](int i, int j, int k) { return i; };

    // Completion well (ISourceSink Abstraction)
    std::vector<std::shared_ptr<ISourceSink>> sources;
    sources.push_back(std::make_shared<ConstantRateWell>(nx / 2, 0, 0, 0, 200.0, q_scale, idx_func));

    Spatial1D spatial(nx, dx);
    auto state = std::make_shared<Reservoir1DState>(spatial, initial_p);
    auto model = std::make_shared<Reservoir1DModel>(k, phi, mu, ct, B, area, sources);

    // Create a well location mask for visualization
    std::vector<double> well_mask(nx, 0.0);
    well_mask[well_idx] = 1.0;
    
    // Using Tridiagonal solver for 1D for efficiency and stability
    auto solver = std::make_shared<LinearTridiagonalSolver>();
    auto integrator = std::make_shared<ImplicitEulerIntegrator>();

    StandardSimulator sim(model, state, solver, integrator);

    std::cout << "Starting 1D Reservoir depletion simulation...\n";
    if (enable_vtk) std::cout << "VTK Export: Enabled\n";
    if (enable_csv) std::cout << "CSV Export: Enabled\n";

    std::ofstream csv;
    if (enable_csv) {
        csv.open("exports/reservoir_results.csv");
        csv << "time,index,x,pressure\n";
    }

    auto logger = [&csv, nx, dx, well_idx, &well_mask, enable_vtk, enable_csv](double t, const IState& s) {
        static int step = 0;
        const auto& r_state = dynamic_cast<const Reservoir1DState&>(s);
        
        if (step % 20 == 0) {
            std::cout << "Time: " << std::fixed << std::setprecision(1) << t 
                      << " hr | Well Pressure: " << r_state.pressures[well_idx] << " psi\n";
        }
        
        if (step % 10 == 0) {
            if (enable_vtk) {
                char filename[100];
                std::sprintf(filename, "exports/reservoir_1d_%03d.vti", step / 10);
                
                std::vector<VTKField> fields = {
                    {"Pressure", r_state.pressures},
                    {"WellLocation", well_mask}
                };
                VTKExporter::export_vti_multi_1d(filename, nx, dx, fields);
            }

            if (enable_csv && csv.is_open()) {
                for (int i = 0; i < r_state.spatial.nx; ++i) {
                    csv << t << "," << i << "," << (i + 0.5) * dx << "," << r_state.pressures[i] << "\n";
                }
            }
        }
        step++;
    };

    sim.run(t_end, dt, logger);

    if (enable_csv && csv.is_open()) csv.close();
    std::cout << "Simulation Successful.\n";

    return 0;
}
