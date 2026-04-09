#include <iostream>
#include <vector>
#include <memory>
#include <iomanip>
#include "state.hpp"
#include "model.hpp"
#include "lib/solvers.hpp"
#include "lib/integrators.hpp"

using namespace numerical_methods;
using namespace numerical_methods::physics_heat;

int main() {
    int nx = 61, ny = 61;
    double dx = 0.1, dy = 0.1, alpha = 0.1, dt = 0.1, t_end = 2.0;

    auto state = std::make_shared<Heat2DState>(nx, ny, dx, dy, 20.0);
    for (int j = 20; j <= 40; ++j) {
        for (int i = 20; i <= 40; ++i) {
            state->temperatures[state->idx(i, j)] = 100.0;
        }
    }

    auto model = std::make_shared<Heat2DModel>(alpha, 20.0, 20.0, 20.0, 20.0);
    auto solver = std::make_shared<ConjugateGradientSolver>();
    auto integrator = std::make_shared<ImplicitEulerIntegrator>();

    StandardSimulator sim(model, state, solver, integrator);

    std::cout << "Starting 2D Implicit Heat Simulation (Generic Simulator)\n";

    auto grid_logger = [](double t, const IState& s) {
        static int step_count = 0;
        if (step_count % 5 == 0) {
            const auto& h_state = dynamic_cast<const Heat2DState&>(s);
            std::cout << "Time: " << std::fixed << std::setprecision(2) << t << "\n";
            
            // Print 10x10 subsample
            int rows = 10, cols = 10;
            int step_x = h_state.nx / cols;
            int step_y = h_state.ny / rows;

            for (int j = 0; j < h_state.ny; j += step_y) {
                for (int i = 0; i < h_state.nx; i += step_x) {
                    double T = h_state.temperatures[h_state.idx(i, j)];
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
