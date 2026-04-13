#include <omp.h>
#include <iostream>
#include <vector>
#include <memory>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "state.hpp"
#include "model.hpp"
#include "lib/integrators.hpp"
#include "lib/modules.hpp"
#include "lib/config_reader.hpp"

using namespace num;
using namespace mod;
using namespace top;
using namespace mod::wave;

int main(int argc, char** argv) {
    std::string config_file = "input/wave_1d.txt";
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

    int nx = config.get("nx", 100);
    double dx = config.get("dx", 0.1);
    double c = config.get("c", 1.0);
    double t_end = config.get("t_end", 10.0);
    double dt = config.get("dt", 0.01);

    Spatial1D spatial(nx, dx);
    auto state = std::make_shared<Wave1DState>(spatial, 0.0);
    
    // Initial guess: Gaussian pulse in the center
    for (int i = 0; i < nx; ++i) {
        double x = i * spatial.dx;
        state->u[i] = std::exp(-std::pow(x - 5.0, 2) / 0.5);
    }

    auto model = std::make_shared<Wave1DModel>(c);
    auto integrator = std::make_shared<RungeKutta4Integrator>();

    StandardSimulator sim(model, state, nullptr, integrator);

    std::cout << "Starting 1D Wave Simulation\n";
    if (enable_csv) std::cout << "CSV Export: Enabled\n";

    std::ofstream csv;
    if (enable_csv) {
        csv.open("exports/wave_results.csv");
        csv << "time,index,x,u,v\n";
    }

    auto logger = [&csv, enable_csv](double t, const IState& s) {
        static int step = 0;
        const auto& w_state = dynamic_cast<const Wave1DState&>(s);
        
        if (step % 50 == 0) {
            std::cout << "Time: " << std::fixed << std::setprecision(3) << t 
                      << " | Mid-point Displacement: " << w_state.u[w_state.u.size() / 2] << "\n";
        }
        
        if (enable_csv && csv.is_open() && step % 10 == 0) {
            for (size_t i = 0; i < w_state.u.size(); ++i) {
                csv << t << "," << i << "," << i * w_state.spatial.dx << "," << w_state.u[i] << "," << w_state.v[i] << "\n";
            }
        }
        step++;
    };

    sim.run(t_end, dt, logger);

    if (enable_csv && csv.is_open()) csv.close();
    std::cout << "Simulation Successful.\n";
    return 0;
}
