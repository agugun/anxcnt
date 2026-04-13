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
    std::string config_file = "input/reservoir_2d.txt";
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

    // Reservoir Parameters
    int nx = config.get("nx", 51);
    int ny = config.get("ny", 51);
    double total_x = config.get("total_x", 2500.0);
    double total_y = config.get("total_y", 2500.0);
    double dx = total_x / nx, dy = total_y / ny;
    double h = config.get("h", 50.0);
    
    double k = config.get("k", 100.0);
    double phi = config.get("phi", 0.2);
    double mu = config.get("mu", 1.0);
    double ct = config.get("ct", 1e-5);
    double B = config.get("B", 1.0);
    double initial_p = config.get("initial_p", 5000.0);
    double dt = config.get("dt", 24.0);
    double t_end = config.get("t_end", 720.0);

    // Producer
    int well_i = nx / 2;
    int well_j = ny / 2;

    // Create a well location mask for visualization
    std::vector<double> well_mask(nx * ny, 0.0);
    well_mask[well_j * nx + well_i] = 1.0;

    // Well term conversion
    double q_scale = (5.615 / 24.0) * B / (phi * ct * h * dx * dy);
    auto idx_func = [nx](int i, int j, int k) { return j * nx + i; };

    std::vector<std::shared_ptr<ISourceSink>> sources;
    sources.push_back(std::make_shared<ConstantRateWell>(well_i, well_j, 0, 0, 500.0, q_scale, idx_func));

    Spatial2D spatial(nx, ny, dx, dy);
    auto state = std::make_shared<Reservoir2DState>(spatial, initial_p);
    auto model = std::make_shared<Reservoir2DModel>(k, phi, mu, ct, B, h, sources);
    
    auto solver = std::make_shared<ConjugateGradientSolver>();
    auto integrator = std::make_shared<ImplicitEulerIntegrator>();

    StandardSimulator sim(model, state, solver, integrator);

    std::cout << "Starting 2D Reservoir Simulation\n";
    if (enable_vtk) std::cout << "VTK Export: Enabled\n";

    auto vtk_logger = [nx, ny, dx, dy, well_i, well_j, &well_mask, enable_vtk](double t, const IState& s) {
        static int step_count = 0;
        const auto& r_state = dynamic_cast<const Reservoir2DState&>(s);
        
        if (enable_vtk && step_count % 5 == 0) {
            char filename[100];
            std::sprintf(filename, "exports/reservoir_2d_%03d.vti", step_count / 5);
            
            std::vector<VTKField> fields = {
                {"Pressure", r_state.pressures},
                {"WellLocation", well_mask}
            };
            VTKExporter::export_vti_multi_2d(filename, nx, ny, dx, dy, fields);
        }

        if (step_count % 5 == 0) {
            std::cout << "Time: " << std::fixed << std::setprecision(1) << t 
                      << " hr | Well Pressure: " << r_state.pressures[r_state.idx(well_i, well_j)] << " psi\n";
        }
        step_count++;
    };

    sim.run(t_end, dt, vtk_logger);

    std::cout << "Simulation Successful.\n";
    return 0;
}
