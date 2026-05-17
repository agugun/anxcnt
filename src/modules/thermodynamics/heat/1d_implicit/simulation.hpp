/**
 * @file simulation.hpp
 * @brief Specialized Simulation class for 1D Implicit Heat Conduction.
 */
#pragma once
#include "lib/simulation.hpp"
#include "lib/linearizers.hpp"
#include "lib/engine_infra.hpp"
#include "lib/discretization.hpp"
#include "lib/solvers.hpp"
#include "lib/integrators.hpp"
#include "state.hpp"
#include "model.hpp"
#include <cmath>

namespace sim {

/**
 * @brief 1D Implicit Heat Simulation.
 */
class HeatSimulation : public Simulation {
public:
    void build(const utl::ConfigReader& config) override {
        size_t nx = config.get("nx", 100);
        double dx = config.get("dx", 0.01);
        double k = config.get("k", config.get("alpha", 0.1));
        double rho = config.get("rho", 1.0);
        double cp = config.get("cp", 1.0);
        double area = config.get("area", 1.0);

        // 1. Grid
        this->grd = std::make_shared<Spatial1D>(nx, dx);

        // 2. Model & Discretizer
        auto cond = num::discretization::heat_cond_1d(nx, dx, k, area);
        Vector storage = num::discretization::heat_storage(nx, dx * area, rho, cp);

        this->mdl = std::make_shared<Heat1DModel>(
            cond, storage, config.get("t_left", 0.0), config.get("t_right", 0.0)
        );
        this->discretizer = std::make_shared<Heat1DDiscretizer>();

        // 3. Numerical Tools
        this->timer = std::make_shared<num::ImplicitEulerIntegrator>();
        this->linearizer = std::make_shared<num::NewtonRaphson>(1e-6, 12, true);
        this->solver = std::make_shared<num::LinearTridiagonalSolver>();
        this->parallel = std::make_shared<utl::SerialParallelManager>();

        if (this->linearizer) this->linearizer->set_sources(this->sources);
    }

    std::unique_ptr<IState> create_initial_state(const utl::ConfigReader& config) {
        size_t nx = config.get("nx", 100);
        double dx = config.get("dx", 0.01);

        auto st = std::make_unique<Heat1DImplicitState>(
            std::static_pointer_cast<Spatial1D>(this->grd),
            0.0
        );

        // Apply Gaussian initial condition
        for (size_t i = 0; i < nx; ++i) {
            double x = i * dx;
            st->temperatures[i] = std::exp(-std::pow(x - 0.5 * (nx-1) * dx, 2) / 0.02);
        }

        return st;
    }
};

} // namespace sim
