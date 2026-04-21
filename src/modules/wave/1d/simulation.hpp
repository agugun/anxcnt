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

namespace mod::wave {

class Wave1DImplicitSimulation {
public:
    struct BuildResult {
        std::unique_ptr<top::SimulationEngine> engine;
        std::unique_ptr<top::IState> initial_state;
        std::shared_ptr<utl::StandardLogger> logger;
    };

    static BuildResult build(const utl::ConfigReader& config) {
        size_t nx = config.get("nx", 200);
        double dx = config.get("dx", 0.01);
        double c = config.get("c", 1.0);
        double rho = config.get("rho", 1.0);
        double area = config.get("area", 1.0);
        
        // 1. Grid and State
        auto spatial = std::make_shared<Spatial1D>(nx, dx);
        auto state = std::make_unique<Wave1DState>(spatial, 0.0);
        
        // Initial condition: Gaussian pulse in displacement
        for (size_t i = 0; i < nx; ++i) {
            double x = i * dx;
            state->u[i] = std::exp(-std::pow(x - 0.5, 2) / 0.01);
        }

        // 2. Physics Model and Discretization
        // In wave eq, Conductance = c^2 * area / dx ? 
        // Actually the discretization utility should handle this.
        auto cond = num::discretization::heat_cond_1d(nx, dx, c * c, area); 
        Vector storage = num::discretization::heat_storage(nx, dx * area, rho, 1.0);
        
        auto model = std::make_shared<Wave1DModel>(cond, storage);
        auto discretizer = std::make_shared<Wave1DDiscretizer>();
        
        // 3. Engine Components
        auto timer = std::make_shared<num::ImplicitEulerIntegrator>();
        auto linearizer = std::make_shared<top::NewtonRaphson>(1e-4, 12, true);
        
        // Use BiCGSTAB as the system is 2N and non-tridiagonal (due to u-v coupling)
        auto solver = std::make_shared<num::BiCGSTABSolver>();
        
        auto pm = std::make_shared<top::SerialParallelManager>();

        auto engine = std::make_unique<top::SimulationEngine>(spatial, model, discretizer, timer, linearizer, solver, pm);

        // 4. Logger / Observer setup
        auto logger = std::make_shared<utl::StandardLogger>(config);
        logger->set_grid(nx, 1, 1, dx);
        logger->add_field("Displacement", [](const top::IState& s) {
            auto v = s.to_vector();
            return Vector(v.begin(), v.begin() + v.size() / 2);
        });
        
        engine->add_observer(logger);

        return { std::move(engine), std::move(state), logger };
    }
};

} // namespace mod::wave
