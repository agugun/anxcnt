#pragma once
#include "lib/simulation.hpp"
#include "lib/linearizers.hpp"
#include "lib/engine_infra.hpp"
#include "lib/discretization.hpp"
#include "lib/solvers.hpp"
#include "lib/integrators.hpp"
#include "state.hpp"
#include "model.hpp"

namespace sim {

class FluidSimulation : public Simulation {
public:
    void build(const utl::ConfigReader& config) override {
        double L = config.get("L", 1.0);
        double H = config.get("H", 1.0);
        int nx = config.get("nx", 10);
        int ny = config.get("ny", 10);
        double mu = config.get("mu", 0.01);
        double rho = config.get("rho", 1.0);

        // 1. Grid (Mesh)
        auto mesh = std::make_shared<Mesh>(Mesh::generate_quad_mesh(L, H, nx, ny));
        this->grd = mesh;

        // 2. Physics Model and Discretization
        auto mdl_fluid = std::make_shared<FluidModel>(mesh, mu, rho);

        // Setup BCs (Lid-driven cavity example)
        for (int i = 0; i < mesh->num_nodes(); ++i) {
            double x = mesh->nodes[i].x;
            double y = mesh->nodes[i].y;

            // Top lid (y=H)
            if (std::abs(y - H) < 1e-6) {
                mdl_fluid->set_velocity_bc(i, 1.0, 0.0);
            }
            // Walls (x=0, x=L, y=0)
            else if (std::abs(x) < 1e-6 || std::abs(x - L) < 1e-6 || std::abs(y) < 1e-6) {
                mdl_fluid->set_velocity_bc(i, 0.0, 0.0);
            }
        }

        this->mdl = mdl_fluid;
        this->discretizer = std::make_shared<FluidDiscretizer>();

        // 3. Engine Components
        this->timer = std::make_shared<num::ImplicitEulerIntegrator>();
        this->linearizer = std::make_shared<num::NewtonRaphson>(1e-4, 10, true);

        auto bicg = std::make_shared<num::BiCGSTABSolver>();
        bicg->verbose = false;
        this->solver = bicg;

        this->parallel = std::make_shared<utl::SerialParallelManager>();
    }

    std::unique_ptr<IState> create_initial_state(const utl::ConfigReader& config) {
        auto mesh = std::dynamic_pointer_cast<Mesh>(this->grd);
        return std::make_unique<FluidState>(mesh);
    }
};

} // namespace sim
