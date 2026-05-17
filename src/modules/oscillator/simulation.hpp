/**
 * @file simulation.hpp
 * @brief Specialized Simulation class for the Harmonic Oscillator.
 */
#pragma once
#include "lib/simulation.hpp"
#include "lib/integrators.hpp"
#include "lib/linearizers.hpp"
#include "lib/solvers.hpp"
#include "lib/engine_infra.hpp"
#include "state.hpp"
#include "model.hpp"

namespace sim {

/**
 * @brief Custom Integrator to enforce a specific dt.
 */
class FixedDtIntegrator : public num::ImplicitEulerIntegrator {
    double fixed_dt;
public:
    FixedDtIntegrator(double dt) : fixed_dt(dt) {}
    double compute_dt(const IState& st, double t) const override {
        return fixed_dt;
    }
};

/**
 * @brief Specialized Simulation for the Harmonic Oscillator.
 */
class OscillatorSimulation : public Simulation {
public:
    void build(const utl::ConfigReader& config) override {
        double m = config.get("m", 1.0);
        double c = config.get("c", 0.0);
        double k = config.get("k", 4.0);
        double dt = config.get("dt", 0.05);

        // 1. Core Physics
        this->grd = std::make_shared<OscillatorGrid>();
        this->mdl = std::make_shared<OscillatorModel>(m, c, k);
        this->discretizer = std::make_shared<OscillatorDiscretizer>();

        // 2. Numerical Components
        this->timer = std::make_shared<FixedDtIntegrator>(dt);
        this->linearizer = std::make_shared<num::NewtonRaphson>(1e-6, 10, false);
        this->solver = std::make_shared<num::LUSolver>();
        this->parallel = std::make_shared<utl::SerialParallelManager>();

        if (this->linearizer) this->linearizer->set_sources(this->sources);
    }

    std::unique_ptr<IState> create_initial_state(const utl::ConfigReader& config) {
        return std::make_unique<OscillatorState>(
            config.get("x0", 1.0),
            config.get("v0", 0.0)
        );
    }
};

} // namespace sim
