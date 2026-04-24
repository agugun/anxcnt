#pragma once
#include "lib/simulation_engine.hpp"
#include "lib/linearizers.hpp"
#include "lib/engine_infra.hpp"
#include "lib/discretization.hpp"
#include "lib/utils/config_reader.hpp"
#include "lib/utils/logger.hpp"
#include "state.hpp"
#include "model.hpp"
#include "discretizer.hpp"
#include "lib/solvers.hpp"
#include "lib/integrators.hpp"

namespace mod::pressure {

class Pressure1DImplicitSimulation {
public:
    struct BuildResult {
        std::unique_ptr<top::SimulationEngine> engine;
        std::unique_ptr<top::IState> initial_state;
        std::shared_ptr<utl::StandardLogger> logger;
    };

    static BuildResult build(const utl::ConfigReader& config) {
        size_t nx = config.get("nx", 100);
        double dx = config.get("dx", 10.0);
        double k = config.get("k", 50.0); // mD
        double phi = config.get("phi", 0.2);
        double mu = config.get("mu", 1.0); // cP
        double ct = config.get("ct", 1e-6); // 1/psi
        double area = config.get("area", 100.0);
        
        // 1. Grid and State
        auto spatial = std::make_shared<Spatial1D>(nx, dx);
        auto state = std::make_unique<Pressure1DState>(spatial, config.get("p_initial", 2000.0));
        
        // 2. Physics Model and Discretization
        auto cond = num::discretization::pressure_cond_1d(nx, dx, k, mu, area);
        Vector storage = num::discretization::pressure_storage(nx, dx * area, phi, ct);
        
        auto model = std::make_shared<Pressure1DModel>(
            cond, storage, config.get("p_left", 3000.0), config.get("p_right", 1000.0)
        );
        auto discretizer = std::make_shared<Pressure1DDiscretizer>();
        
        // 3. Engine Components
        auto timer = std::make_shared<num::ImplicitEulerIntegrator>();
        auto linearizer = std::make_shared<num::NewtonRaphson>(1e-6, 12, true);
        auto solver = std::make_shared<num::LinearTridiagonalSolver>();
        auto pm = std::make_shared<top::SerialParallelManager>();

        auto engine = std::make_unique<top::SimulationEngine>(spatial, model, discretizer, timer, linearizer, solver, pm);

        // 4. Logger / Observer setup
        auto logger = std::make_shared<utl::StandardLogger>(config);
        logger->set_grid(nx, 1, 1, dx);
        logger->add_field("Pressure", [](const top::IState& s) {
            return s.to_vector();
        });
        
        engine->add_observer(logger);

        return { std::move(engine), std::move(state), logger };
    }
};

} // namespace mod::pressure
