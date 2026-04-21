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

namespace mod::fluid {

class FluidSimulationBuilder {
public:
    struct BuildResult {
        std::unique_ptr<top::SimulationEngine> engine;
        std::unique_ptr<top::IState> initial_state;
        std::shared_ptr<utl::StandardLogger> logger;
    };

    static BuildResult build(const utl::ConfigReader& config) {
        double L = config.get("L", 1.0);
        double H = config.get("H", 1.0);
        int nx = config.get("nx", 10);
        int ny = config.get("ny", 10);
        double mu = config.get("mu", 0.01);
        double rho = config.get("rho", 1.0);
        
        // 1. Grid (Mesh) and State
        auto mesh = std::make_shared<Mesh>(Mesh::generate_quad_mesh(L, H, nx, ny));
        auto state = std::make_unique<FluidState>(mesh);
        
        // 2. Physics Model and Discretization
        auto model = std::make_shared<FluidModel>(mesh, mu, rho);
        auto discretizer = std::make_shared<FluidDiscretizer>();
        
        // Setup BCs (Lid-driven cavity example)
        for (int i = 0; i < mesh->num_nodes(); ++i) {
            double x = mesh->nodes[i].x;
            double y = mesh->nodes[i].y;
            
            // Top lid (y=H)
            if (std::abs(y - H) < 1e-6) {
                model->set_velocity_bc(i, 1.0, 0.0);
            }
            // Walls (x=0, x=L, y=0)
            else if (std::abs(x) < 1e-6 || std::abs(x - L) < 1e-6 || std::abs(y) < 1e-6) {
                model->set_velocity_bc(i, 0.0, 0.0);
            }
        }

        // 3. Engine Components
        auto timer = std::make_shared<num::ImplicitEulerIntegrator>();
        auto linearizer = std::make_shared<top::NewtonRaphson>(1e-4, 10, true);
        
        // FEM systems are usually non-symmetric and larger, BiCGSTAB is standard
        auto solver = std::make_shared<num::BiCGSTABSolver>();
        solver->verbose = false;
        
        auto pm = std::make_shared<top::SerialParallelManager>();

        auto engine = std::make_unique<top::SimulationEngine>(nullptr, model, discretizer, timer, linearizer, solver, pm);

        // 4. Logger / Observer setup
        auto logger = std::make_shared<utl::StandardLogger>(config);
        logger->set_grid(nx, ny, 1, L/nx);
        logger->add_field("U-Velocity", [](const top::IState& s) {
            const auto& fs = static_cast<const FluidState&>(s);
            return fs.u;
        });
        logger->add_field("P-Field", [](const top::IState& s) {
            const auto& fs = static_cast<const FluidState&>(s);
            return fs.p;
        });
        
        engine->add_observer(logger);

        return { std::move(engine), std::move(state), logger };
    }
};

} // namespace mod::fluid
