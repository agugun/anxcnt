#include "lib/simulation_engine.hpp"
#include "lib/linearizers.hpp"
#include "lib/engine_infra.hpp"
#include "lib/discretization.hpp"
#include "lib/utils/config_reader.hpp"
#include "lib/utils/logger.hpp"
#include "state.hpp"
#include "model.hpp"
#include "lib/integrators.hpp"
#include "lib/solvers.hpp"

namespace mod::heat {

class Heat1DImplicitSimulation {
public:
    struct BuildResult {
        std::unique_ptr<top::SimulationEngine> engine;
        std::unique_ptr<top::IState> initial_state;
        std::shared_ptr<utl::StandardLogger> logger;
    };

    static BuildResult build(const utl::ConfigReader& config) {
        size_t nx = config.get("nx", 100);
        double dx = config.get("dx", 0.01);
        double k = config.get("k", config.get("alpha", 0.1));
        double rho = config.get("rho", 1.0);
        double cp = config.get("cp", 1.0);
        double area = config.get("area", 1.0);
        
        // 1. Grid and State
        auto spatial = std::make_shared<Spatial1D>(nx, dx);
        auto state = std::make_unique<Heat1DImplicitState>(spatial, 0.0);
        for (size_t i = 0; i < nx; ++i) {
            double x = i * dx;
            state->temperatures[i] = std::exp(-std::pow(x - 0.5 * (nx-1) * dx, 2) / 0.02);
        }

        // 2. Physics Model and Discretization
        auto cond = num::discretization::heat_cond_1d(nx, dx, k, area);
        Vector storage = num::discretization::heat_storage(nx, dx * area, rho, cp);
        auto model = std::make_shared<Heat1DModel>(cond, storage, config.get("t_left", 0.0), config.get("t_right", 0.0));
        auto discretizer = std::make_shared<Heat1DDiscretizer>();
        
        // 3. Engine Components
        auto timer = std::make_shared<num::ImplicitEulerIntegrator>();
        auto linearizer = std::make_shared<top::NewtonRaphson>(1e-6, 12, true);
        auto solver = std::make_shared<num::LinearTridiagonalSolver>();
        auto pm = std::make_shared<top::SerialParallelManager>();

        auto engine = std::make_unique<top::SimulationEngine>(spatial, model, discretizer, timer, linearizer, solver, pm);

        // 4. Logger / Observer setup
        auto logger = std::make_shared<utl::StandardLogger>(config);
        logger->set_grid(nx, 1, 1, dx);
        logger->add_field("Temperature", [](const top::IState& s) {
            return s.to_vector();
        });
        
        engine->add_observer(logger);

        return { std::move(engine), std::move(state), logger };
    }
};

} // namespace mod::heat
