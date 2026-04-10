#include <iostream>
#include <vector>
#include <memory>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "state.hpp"
#include "model.hpp"
#include "lib/integrators.hpp"

using namespace num;
using namespace mod;
using namespace top;
using namespace mod::physics_wave;

int main() {
    int nx = 201;
    double dx = 0.01;
    double c = 1.0;
    double dt = 0.005;
    double t_end = 2.0;

    auto state = std::make_shared<Wave1DState>(nx, dx, 0.0);
    for (int i = 0; i < nx; ++i) {
        state->u[i] = std::exp(-std::pow(i * dx - 1.0, 2) / 0.02);
    }

    auto model = std::make_shared<Wave1DModel>(c);
    auto integrator = std::make_shared<RungeKutta4Integrator>();

    StandardSimulator sim(model, state, nullptr, integrator);

    std::cout << "Starting 1D Wave Simulation (Generic Simulator)\n";
    std::cout << "Results saved to: exports/wave_results.csv\n";

    std::ofstream csv("exports/wave_results.csv");
    csv << "time,index,x,u,v\n";

    auto logger = [&csv](double t, const IState& s) {
        static int step = 0;
        const auto& w_state = dynamic_cast<const Wave1DState&>(s);
        
        if (step % 50 == 0) {
            std::cout << "Time: " << std::fixed << std::setprecision(3) << t 
                      << " | Mid-point Displacement: " << w_state.u[w_state.u.size() / 2] << "\n";
        }
        
        if (step % 10 == 0) {
            for (size_t i = 0; i < w_state.u.size(); ++i) {
                csv << t << "," << i << "," << i * w_state.dx << "," << w_state.u[i] << "," << w_state.v[i] << "\n";
            }
        }
        step++;
    };

    sim.run(t_end, dt, logger);

    csv.close();
    std::cout << "Simulation Successful.\n";
    return 0;
}
