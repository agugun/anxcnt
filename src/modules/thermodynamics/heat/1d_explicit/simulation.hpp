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

class Heat1DExplicitSimulation : public Simulation {
public:
    void build(const utl::ConfigReader& config) override {
        size_t nx = config.get("nx", 100);
        double dx = config.get("dx", 0.01);
        double k = config.get("k", 0.1);
        double rho = config.get("rho", 1.0);
        double cp = config.get("cp", 1.0);
        double area = config.get("area", 1.0);

        // 1. Grid
        this->grd = std::make_shared<Spatial1D>(nx, dx);

        // 2. Physics Model and Discretization
        auto cond = num::discretization::heat_cond_1d(nx, dx, k, area);
        std::vector<double> storage = num::discretization::heat_storage(nx, dx * area, rho, cp);

        this->mdl = std::make_shared<Heat1DExplicitModel>(cond, storage);
        this->discretizer = std::make_shared<Heat1DExplicitDiscretizer>();

        // 3. Engine Components (Explicit)
        this->timer = std::make_shared<num::ForwardEulerIntegrator>();
        this->linearizer = std::make_shared<num::ExplicitLinearizer>();
        this->solver = std::make_shared<num::LinearTridiagonalSolver>();
        this->parallel = std::make_shared<utl::SerialParallelManager>();
    }

    std::unique_ptr<IState> create_initial_state(const utl::ConfigReader& config) {
        size_t nx = config.get("nx", 100);
        double dx = config.get("dx", 0.01);
        auto spatial = std::dynamic_pointer_cast<Spatial1D>(this->grd);
        auto st = std::make_unique<Heat1DExplicitState>(*spatial, config.get("initial_temp", 25.0));
        st->temperatures[nx / 2] = 100.0; // Initial condition splash
        return st;
    }
};

} // namespace sim
