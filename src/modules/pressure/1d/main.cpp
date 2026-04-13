#include <iostream>
#include <vector>
#include <memory>
#include <fstream>
#include <iomanip>
#include "state.hpp"
#include "model.hpp"
#include "lib/solvers.hpp"
#include "lib/integrators.hpp"
#include "lib/modules.hpp"
#include "lib/config_reader.hpp"

using namespace num;
using namespace mod;
using namespace top;
using namespace mod::pressure;

int main(int argc, char** argv) {
    std::string config_file = "input/pressure_1d.txt";
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] != '-') {
            config_file = argv[i];
            break;
        }
    }

    ConfigReader config;
    config.load(config_file);

    bool enable_csv = config.get("enable_csv", 0) != 0;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--csv") enable_csv = true;
    }

    int nx = config.get("nx", 50);
    double dx = config.get("dx", 20.0);
    double k = config.get("k", 50.0);
    double phi = config.get("phi", 0.25);
    double mu = config.get("mu", 2.0);
    double ct = config.get("ct", 1e-5);
    double p_left = config.get("p_left", 1500.0);
    double p_right = config.get("p_right", 3000.0);

    Spatial1D spatial(nx, dx);
    auto state = std::make_shared<Pressure1DState>(spatial, p_right);
    state->pressures[0] = p_left;
    state->pressures[nx-1] = p_right;

    auto model = std::make_shared<Pressure1DModel>(k, phi, mu, ct, p_left, p_right);
    auto solver = std::make_shared<LinearTridiagonalSolver>();
    auto integrator = std::make_shared<ImplicitEulerIntegrator>();

    StandardSimulator sim(model, state, solver, integrator);

    std::cout << "Starting 1D Pressure Diffusivity Simulation\n";
    if (enable_csv) std::cout << "CSV Export: Enabled\n";

    std::ofstream csv;
    if (enable_csv) {
        csv.open("exports/pressure_results.csv");
        csv << "time,index,x,pressure\n";
    }

    auto logger = [&csv, enable_csv](double t, const IState& s) {
        static int step = 0;
        const auto& p_state = dynamic_cast<const Pressure1DState&>(s);
        
        if (step % 20 == 0) {
            std::cout << "Time: " << std::fixed << std::setprecision(1) << t 
                      << " hr | Mid-point Pressure: " << p_state.pressures[p_state.pressures.size() / 2] << " psi\n";
        }
        
        if (enable_csv && csv.is_open() && step % 5 == 0) {
            for (size_t i = 0; i < p_state.pressures.size(); ++i) {
                csv << t << "," << i << "," << i * p_state.spatial.dx << "," << p_state.pressures[i] << "\n";
            }
        }
        step++;
    };

    sim.run(2400.0, 24.0, logger);

    if (enable_csv && csv.is_open()) csv.close();
    std::cout << "Simulation Successful.\n";
    return 0;
}
