#include <iostream>
#include <vector>
#include <memory>
#include <fstream>
#include <iomanip>
#include "state.hpp"
#include "model.hpp"
#include "lib/integrators.hpp"

using namespace num;
using namespace mod;
using namespace top;
using namespace mod::physics_mba;

int main() {
    double pi = 5000.0, N = 100e6, ct = 1e-5, q = 5000.0;
    
    auto state = std::make_shared<MBState>(pi);
    auto model = std::make_shared<MBModel>(N, ct, q);
    auto integrator = std::make_shared<ForwardEulerIntegrator>();

    StandardSimulator sim(model, state, nullptr, integrator);

    std::cout << "Starting Material Balance (MBA) Simulation (Generic Simulator)\n";

    std::ofstream csv("exports/mba_results.csv");
    csv << "time,pressure,cum_production\n";

    auto logger = [&csv, pi](double t, const IState& s) {
        static int step = 0;
        const auto& mb_state = dynamic_cast<const MBState&>(s);
        
        if (step % 50 == 0) {
            std::cout << "Day: " << std::setw(3) << (int)t 
                      << " | Pressure: " << std::fixed << std::setprecision(1) 
                      << mb_state.pressure << " psi\n";
        }
        
        if (step % 10 == 0) {
            double cum_prod = (pi - mb_state.pressure) * 1e6; // Dummy scaling
            csv << t << "," << mb_state.pressure << "," << cum_prod << "\n";
        }
        step++;
    };

    sim.run(365.0, 1.0, logger);

    csv.close();
    std::cout << "Simulation Successful.\n";
    return 0;
}
