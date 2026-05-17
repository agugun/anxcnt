/**
 * @file simulation.hpp
 * @brief Specialized Simulation class for 1D Pressure Diffusivity.
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

namespace sim {

/**
 * @brief 1D Pressure Simulation.
 * Implements the build factory method to assemble a pressure diffusivity solver.
 */
class PressureSimulation : public Simulation {
public:
    void build(const utl::ConfigReader& config) override {
        size_t nx = config.get("nx", 100);
        double dx = config.get("dx", 10.0);
        double k = config.get("k", 50.0);
        double phi = config.get("phi", 0.2);
        double mu = config.get("mu", 1.0);
        double ct = config.get("ct", 1e-6);
        double area = config.get("area", 100.0);

        // 1. Grid
        this->grd = std::make_shared<geo::Spatial1D>(nx, dx);

        // 2. Model & Discretizer
        auto cond = disc::pressure_cond_1d(nx, dx, k, mu, area);
        Vector storage = disc::pressure_storage(nx, dx * area, phi, ct);

        this->mdl = std::make_shared<Pressure1DModel>(
            cond, storage, config.get("p_left", 3000.0), config.get("p_right", 1000.0)
        );
        this->discretizer = std::make_shared<Pressure1DDiscretizer>();

        // 3. Numerical Tools
        this->timer = std::make_shared<num::ImplicitEulerIntegrator>();
        this->linearizer = std::make_shared<num::NewtonRaphson>(1e-6, 12, true);
        this->solver = std::make_shared<num::LinearTridiagonalSolver>();
        this->parallel = std::make_shared<utl::SerialParallelManager>();

        if (this->linearizer) this->linearizer->set_sources(this->sources);
    }

    /**
     * @brief Helper to create the initial state based on config.
     */
    std::unique_ptr<IState> create_initial_state(const utl::ConfigReader& config) {
        return std::make_unique<Pressure1DState>(
            std::static_pointer_cast<geo::Spatial1D>(this->grd),
            config.get("p_initial", 2000.0)
        );
    }
};

} // namespace sim
