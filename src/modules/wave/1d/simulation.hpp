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

class Wave1DSimulation : public Simulation {
public:
    void build(const utl::ConfigReader& config) override {
        size_t nx = config.get("nx", 200);
        double dx = config.get("dx", 0.01);
        double c = config.get("c", 1.0);
        double rho = config.get("rho", 1.0);
        double area = config.get("area", 1.0);

        // 1. Grid
        this->grd = std::make_shared<Spatial1D>(nx, dx);

        // 2. Physics Model and Discretization
        auto cond = num::discretization::heat_cond_1d(nx, dx, c * c, area);
        Vector storage = num::discretization::heat_storage(nx, dx * area, rho, 1.0);

        this->mdl = std::make_shared<Wave1DModel>(cond, storage);
        this->discretizer = std::make_shared<Wave1DDiscretizer>();

        // 3. Engine Components
        this->timer = std::make_shared<num::ImplicitEulerIntegrator>();
        this->linearizer = std::make_shared<num::NewtonRaphson>(1e-4, 12, true);
        this->solver = std::make_shared<num::BiCGSTABSolver>();
        this->parallel = std::make_shared<utl::SerialParallelManager>();
    }

    std::unique_ptr<IState> create_initial_state(const utl::ConfigReader& config) {
        size_t nx = config.get("nx", 200);
        double dx = config.get("dx", 0.01);
        auto spatial = std::dynamic_pointer_cast<Spatial1D>(this->grd);
        auto st = std::make_unique<Wave1DState>(spatial, 0.0);

        for (size_t i = 0; i < nx; ++i) {
            double x = i * dx;
            st->u[i] = std::exp(-std::pow(x - 0.5, 2) / 0.01);
        }
        return st;
    }
};

} // namespace sim
