#pragma once
#include "lib/simulation_engine.hpp"
#include "lib/linearizers.hpp"
#include "lib/engine_infra.hpp"
#include "lib/discretization.hpp"
#include "lib/utils/config_reader.hpp"
#include "lib/utils/logger.hpp"
#include "state.hpp"
#include "model.hpp"
#include "lib/solvers.hpp"
#include "lib/integrators.hpp"

namespace mod::burgers {

class BurgersSimulationBuilder {
public:
    struct BuildResult {
        std::unique_ptr<top::SimulationEngine> engine;
        std::unique_ptr<top::IState> initial_state;
        std::shared_ptr<utl::StandardLogger> logger;
    };

    static BuildResult build(const utl::ConfigReader& config) {
        size_t nx = config.get("nx", 100);
        double dx = config.get("dx", 0.01);
        double nu = config.get("nu", 0.01);
        
        // 1. Grid and State
        auto spatial = std::make_shared<Spatial1D>(nx, dx);
        auto state = std::make_unique<BurgersState>(spatial, 0.0);
        
        // Initial condition: Sine wave
        for (size_t i = 0; i < nx; ++i) {
            double x = i * dx;
            state->u[i] = std::sin(2.0 * M_PI * x);
        }

        // 2. Physics Model and Discretization
        auto model = std::make_shared<BurgersModel>(nu, dx);
        auto discretizer = std::make_shared<BurgersDiscretizer>();
        
        // 3. Engine Components
        auto timer = std::make_shared<num::ImplicitEulerIntegrator>();
        auto linearizer = std::make_shared<num::NewtonRaphson>(1e-6, 15, true);
        
        // 1D non-linear system, use Tridiagonal if possible or BiCGSTAB
        auto solver = std::make_shared<num::BiCGSTABSolver>();
        solver->verbose = false;
        
        auto pm = std::make_shared<top::SerialParallelManager>();

        auto engine = std::make_unique<top::SimulationEngine>(spatial, model, discretizer, timer, linearizer, solver, pm);

        // 4. Logger / Observer setup
        auto logger = std::make_shared<utl::StandardLogger>(config);
        logger->set_grid(nx, 1, 1, dx);
        logger->add_field("Velocity", [](const top::IState& s) {
            return s.to_vector();
        });
        
        engine->add_observer(logger);

        return { std::move(engine), std::move(state), logger };
    }
};

} // namespace mod::burgers
