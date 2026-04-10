#include <iostream>
#include <vector>
#include <memory>
#include <fstream>
#include <iomanip>
#include "state.hpp"
#include "model.hpp"
#include "lib/solvers.hpp"
#include "lib/integrators.hpp"

using namespace num;
using namespace mod;
using namespace top;
using namespace mod::physics_pressure;

int main() {
    int nx = 50;
    double dx = 20.0;
    double k = 50.0, phi = 0.25, mu = 2.0, ct = 1e-5;
    double p_left = 1500.0, p_right = 3000.0;

    auto state = std::make_shared<Pressure1DState>(nx, dx, p_right);
    state->pressures[0] = p_left;
    state->pressures[nx-1] = p_right;

    auto model = std::make_shared<Pressure1DModel>(k, phi, mu, ct, p_left, p_right);
    auto solver = std::make_shared<LinearTridiagonalSolver>();
    auto integrator = std::make_shared<ImplicitEulerIntegrator>();

    StandardSimulator sim(model, state, solver, integrator);

    std::cout << "Starting 1D Pressure Diffusivity Simulation (Generic Simulator)\n";
    std::cout << "Results saved to: exports/pressure_results.csv\n";

    std::ofstream csv("exports/pressure_results.csv");
    csv << "time,index,x,pressure\n";

    auto logger = [&csv](double t, const IState& s) {
        static int step = 0;
        const auto& p_state = dynamic_cast<const Pressure1DState&>(s);
        
        if (step % 20 == 0) {
            std::cout << "Time: " << std::fixed << std::setprecision(1) << t 
                      << " hr | Mid-point Pressure: " << p_state.pressures[p_state.pressures.size() / 2] << " psi\n";
        }
        
        if (step % 5 == 0) {
            for (size_t i = 0; i < p_state.pressures.size(); ++i) {
                csv << t << "," << i << "," << i * p_state.dx << "," << p_state.pressures[i] << "\n";
            }
        }
        step++;
    };

    sim.run(2400.0, 24.0, logger);

    csv.close();
    std::cout << "Simulation Successful.\n";
    return 0;
}
