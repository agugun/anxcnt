#include <iostream>
#include <fstream>
#include <memory>
#include "state.hpp"
#include "model.hpp"
#include "simulation.hpp"
#include "lib/simulation_engine.hpp"
#include "lib/integrators.hpp"
#include "lib/linearizers.hpp"
#include "lib/solvers.hpp"
#include "lib/engine_infra.hpp"

using namespace mod::oscillator;
using namespace top;
using namespace num;

/**
 * @brief Custom Integrator to enforce a specific dt.
 */
class FixedDtIntegrator : public num::ImplicitEulerIntegrator {
    double fixed_dt;
public:
    FixedDtIntegrator(double dt) : fixed_dt(dt) {}
    double get_next_timestep(const top::IState& state, double t) const override { 
        return fixed_dt; 
    }
};

/**
 * @brief Observer to write state variables to CSV.
 */
class CSVObserver : public IObserver {
    std::ofstream out;
public:
    CSVObserver(const std::string& filename) {
        out.open(filename);
        out << "Time,Position,Velocity\n";
    }
    void on_step_complete(double t, int step, const IState& base_state) override {
        const auto& state = dynamic_cast<const OscillatorState&>(base_state);
        out << t << "," << state.x << "," << state.v << "\n";
    }
};

int main() {
    std::cout << "Starting Harmonic Oscillator Simulation..." << std::endl;

    // 1. Initial State: x(0) = 1, v(0) = 0
    auto initial_state = std::make_unique<OscillatorState>(1.0, 0.0);

    // 2. Model: m=1.0, c=0.0 (undamped to match SymPy phase 1), k=4.0 (omega=2.0)
    auto model = std::make_shared<OscillatorModel>(1.0, 0.0, 4.0);
    auto grid = std::make_shared<OscillatorGrid>();
    auto discretizer = std::make_shared<OscillatorDiscretizer>();

    // 3. Numerical Engine Components
    double dt = 0.05;
    auto timer = std::make_shared<FixedDtIntegrator>(dt);
    auto linearizer = std::make_shared<NewtonRaphson>(1e-6, 10, false);
    auto solver = std::make_shared<LUSolver>();
    auto pm = std::make_shared<SerialParallelManager>();

    SimulationEngine engine(grid, model, discretizer, timer, linearizer, solver, pm);

    auto observer = std::make_shared<CSVObserver>("exports/oscillator_output.csv");
    engine.add_observer(observer);

    // 4. Run Simulation
    double t_max = 10.0;
    
    // We add an epsilon to ensure the final step completes
    engine.simulate(t_max + 1e-4, dt, std::move(initial_state));

    std::cout << "Simulation complete. Output saved to exports/oscillator_output.csv" << std::endl;
    return 0;
}
