#include <omp.h>
#include <iostream>
#include <vector>
#include <memory>
#include <fstream>
#include <iomanip>
#include "state.hpp"
#include "model.hpp"
#include "lib/integrators.hpp"
#include "lib/modules.hpp"
#include "lib/config_reader.hpp"

using namespace num;
using namespace mod;
using namespace top;
using namespace mod::physics_mba;

int main(int argc, char** argv) {
    std::string config_file = "input/reservoir_mba.txt";
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] != '-') {
            config_file = argv[i];
            break;
        }
    }

    ConfigReader config;
    config.load(config_file);
    int num_threads = config.get("num_threads", 1);
    omp_set_num_threads(num_threads);

    bool enable_csv = config.get("enable_csv", 0) != 0;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--csv") enable_csv = true;
    }

    double pi = config.get("pi", 5000.0);
    double N = config.get("N", 100e6);
    double ct = config.get("ct", 1e-5);
    double q = config.get("q", 5000.0);
    
    auto state = std::make_shared<MBState>(pi);
    auto model = std::make_shared<MBModel>(N, ct, q);
    auto integrator = std::make_shared<ForwardEulerIntegrator>();

    StandardSimulator sim(model, state, nullptr, integrator);

    std::cout << "Starting Material Balance (MBA) Simulation\n";
    if (enable_csv) std::cout << "CSV Export: Enabled\n";

    std::ofstream csv;
    if (enable_csv) {
        csv.open("exports/mba_results.csv");
        csv << "time,pressure,cum_production\n";
    }

    auto logger = [&csv, pi, enable_csv](double t, const IState& s) {
        static int step = 0;
        const auto& mb_state = dynamic_cast<const MBState&>(s);
        
        if (step % 50 == 0) {
            std::cout << "Day: " << std::setw(3) << (int)t 
                      << " | Pressure: " << std::fixed << std::setprecision(1) 
                      << mb_state.pressure << " psi\n";
        }
        
        if (enable_csv && csv.is_open() && step % 10 == 0) {
            double cum_prod = (pi - mb_state.pressure) * 1e6; // Dummy scaling
            csv << t << "," << mb_state.pressure << "," << cum_prod << "\n";
        }
        step++;
    };

    double t_end = config.get("t_end", 365.0);
    double dt = config.get("dt", 1.0);
    sim.run(t_end, dt, logger);

    if (enable_csv && csv.is_open()) csv.close();
    std::cout << "Simulation Successful.\n";
    return 0;
}
