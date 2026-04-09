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

using namespace numerical_methods;
using namespace numerical_methods::physics_heat;

int main() {
    int nx = 101;
    double dx = 0.01;
    double alpha = 0.1;
    double dt = 0.01;
    double t_end = 1.0;

    auto state = std::make_shared<Heat1DState>(nx, dx, 0.0);
    for (int i = 0; i < nx; ++i) {
        state->temperatures[i] = std::exp(-std::pow(i * dx - 0.5, 2) / 0.02);
    }

    auto model = std::make_shared<Heat1DModel>(alpha, 0.0, 0.0);
    auto solver = std::make_shared<LinearTridiagonalSolver>();
    auto integrator = std::make_shared<ImplicitEulerIntegrator>();

    StandardSimulator sim(model, state, solver, integrator);

    std::cout << "Starting Implicit Heat Simulation (Generic Simulator)\n";

    std::ofstream csv("exports/heat_results.csv");
    csv << "time,index,x,temp\n";

    auto logger = [&csv](double t, const IState& s) {
        static int step = 0;
        const auto& h_state = dynamic_cast<const Heat1DState&>(s);
        
        // Log to console every 20 steps
        if (step % 20 == 0) {
            std::cout << "Time: " << std::fixed << std::setprecision(3) << t 
                      << " | Center Temp: " << h_state.temperatures[h_state.temperatures.size() / 2] << "\n";
        }
        
        // Log to CSV every 5 steps
        if (step % 5 == 0) {
            for (size_t i = 0; i < h_state.temperatures.size(); ++i) {
                csv << t << "," << i << "," << i * h_state.dx << "," << h_state.temperatures[i] << "\n";
            }
        }
        step++;
    };

    sim.run(t_end, dt, logger);

    csv.close();
    std::cout << "Simulation Successful.\n";
    return 0;
}
