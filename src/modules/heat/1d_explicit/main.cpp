#include <omp.h>
#include <iostream>
#include <vector>
#include <memory>
#include <iomanip>
#include "state.hpp"
#include "model.hpp"
#include "lib/integrators.hpp"
#include "lib/io.hpp"
#include "lib/modules.hpp"
#include "lib/config_reader.hpp"

using namespace num;
using namespace mod;
using namespace top;
using namespace mod::heat;

int main(int argc, char** argv) {
  std::string config_file = "input/heat_1d_explicit.txt";
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

  bool enable_vtk = config.get("enable_vtk", 0) != 0;
  for (int i = 1; i < argc; ++i) {
      if (std::string(argv[i]) == "--vtk") enable_vtk = true;
  }

  size_t num_nodes = config.get("nx", 50);
  double alpha = config.get("alpha", 0.01);
  double dx = config.get("dx", 0.1);
  double dt = config.get("dt", 0.1);
  double t_end = config.get("t_end", 5.0);

  Spatial1D spatial(num_nodes, dx);
  auto state = std::make_shared<Heat1DExplicitState>(spatial, 25.0);
  state->temperatures[num_nodes / 2] = 100.0;

  auto model = std::make_shared<HeatModel>(alpha, spatial);
  auto integrator = std::make_shared<ForwardEulerIntegrator>();

  StandardSimulator sim(model, state, nullptr, integrator);

  std::cout << "Starting Heat Simulation\n";
  if (enable_vtk) std::cout << "VTK Export: Enabled\n";

  auto logger = [num_nodes, dx, enable_vtk](double t, const IState &s) {
    static int step_count = 0;
    const auto &h_state = dynamic_cast<const Heat1DExplicitState &>(s);

    if (enable_vtk) {
        char filename[100];
        std::sprintf(filename, "exports/heat_1d_ex_%03d.vti", step_count);
        VTKExporter::export_vti_1d(filename, h_state.temperatures, (int)h_state.spatial.nx, h_state.spatial.dx, "Temperature");
    }

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