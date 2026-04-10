#include "lib/integrators.hpp"
#include "modules/heat_1d_explicit/model.hpp"
#include "modules/heat_1d_explicit/state.hpp"
#include "lib/io.hpp"
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

using namespace num;
using namespace mod;
using namespace top;
using namespace mod::physics_heat;

int main() {
  size_t num_nodes = 50;
  double alpha = 0.01;
  double dx = 0.1;
  double dt = 0.1;
  double t_end = 5.0;

  auto state = std::make_shared<HeatState>(num_nodes, 25.0);
  state->temperatures[num_nodes / 2] = 100.0;

  auto model = std::make_shared<HeatModel>(alpha, dx);
  auto integrator = std::make_shared<ForwardEulerIntegrator>();

  // Use the generic StandardSimulator
  StandardSimulator sim(model, state, nullptr, integrator);

  std::cout << "Starting Heat Simulation (Callback Logging)\n";

  // Define logging via lambda
  auto logger = [num_nodes, dx](double t, const IState &s) {
    static int step_count = 0;
    const auto &h_state = dynamic_cast<const HeatState &>(s);

    // Export VTI every step
    char filename[100];
    std::sprintf(filename, "exports/heat_1d_ex_%03d.vti", step_count);
    VTKExporter::export_vti_1d(filename, h_state.temperatures, (int)num_nodes, dx, "Temperature");

    if (step_count % 10 == 0) {
      std::cout << "Time: " << std::fixed << std::setprecision(1) << t
                << " | Center Temp: " << std::setprecision(4)
                << h_state.temperatures[h_state.temperatures.size() / 2]
                << "\n";
    }
    step_count++;
  };

  sim.run(t_end, dt, logger);

  std::cout << "Simulation Successful.\n";
  return 0;
}