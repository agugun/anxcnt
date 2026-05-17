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

class BurgersSimulation : public Simulation {
public:
    void build(const utl::ConfigReader& config) override {
        size_t nx = config.get("nx", 100);
        double dx = config.get("dx", 0.01);
        double nu = config.get("nu", 0.01);

        // 1. Grid
        this->grd = std::make_shared<Spatial1D>(nx, dx);

        // 2. Physics Model and Discretization
        this->mdl = std::make_shared<BurgersModel>(nu, dx);
        this->discretizer = std::make_shared<BurgersDiscretizer>();

        // 3. Engine Components
        this->timer = std::make_shared<num::ImplicitEulerIntegrator>();
        this->linearizer = std::make_shared<num::NewtonRaphson>(1e-6, 15, true);

        auto bicg = std::make_shared<num::BiCGSTABSolver>();
        bicg->verbose = false;
        this->solver = bicg;

        this->parallel = std::make_shared<utl::SerialParallelManager>();
    }

    std::unique_ptr<IState> create_initial_state(const utl::ConfigReader& config) {
        size_t nx = config.get("nx", 100);
        double dx = config.get("dx", 0.01);
        auto spatial = std::dynamic_pointer_cast<Spatial1D>(this->grd);
        auto st = std::make_unique<BurgersState>(spatial, 0.0);

        for (size_t i = 0; i < nx; ++i) {
            double x = i * dx;
            st->u[i] = std::sin(2.0 * M_PI * x);
        }
        return st;
    }
};

} // namespace sim
