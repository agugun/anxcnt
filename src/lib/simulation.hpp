/**
 * @file simulation.hpp
 * @brief Template-Based Simulation Framework for Physics Modules.
 */
#pragma once
#include "interfaces.hpp"
#include "utils/config_reader.hpp"
#include <cmath>
#include <vector>
#include <memory>

namespace sim {
using namespace mod;

/**
 * @brief Base class for all physics-specific simulations.
 * Acts as both a Factory (build) and an Orchestrator (run).
 */
class Simulation {
protected:
    std::shared_ptr<geo::IGrid> grd;
    std::shared_ptr<IModel> mdl;
    std::shared_ptr<IDiscretizer> discretizer;
    std::shared_ptr<num::ITimeIntegrator> timer;
    std::shared_ptr<num::ILinearizer> linearizer;
    std::shared_ptr<num::ISolver> solver;
    std::shared_ptr<utl::IParallelManager> parallel;
    std::vector<std::shared_ptr<ISourceSink>> sources;
    std::vector<std::shared_ptr<utl::IObserver>> observers;

public:
    Simulation() = default;
    virtual ~Simulation() = default;

    /**
     * @brief Factory Method: Assembles the module-specific components.
     * Must be implemented by each physics case.
     */
    virtual void build(const utl::ConfigReader& config) = 0;

    void add_observer(std::shared_ptr<utl::IObserver> ob) {
        observers.push_back(ob);
    }

    /**
     * @brief Perform a single time step and notify observers.
     * @return unique_ptr to the next state.
     */
    std::unique_ptr<IState> step(double t, double dt, const IState& st_curr, int step_count = 0) {
        auto st_next = linearizer->resolve(
            st_curr, dt, *grd, *mdl, *discretizer, *timer, *solver, *parallel
        );

        for (auto& obs : observers) {
            obs->on_step_complete(t + dt, step_count, *st_next);
        }

        return st_next;
    }

    /**
     * @brief Run the full simulation until t_max.
     */
    std::unique_ptr<IState> run(double t_max, double dt_initial, std::unique_ptr<IState> st_init) {
        auto st_n = std::move(st_init);
        double t = 0.0;
        int step_count = 0;
        const double time_epsilon = 1e-12 * (std::abs(t_max) + 1.0);

        for (auto& obs : observers) obs->on_simulation_start(*grd);
        for (auto& obs : observers) obs->on_step_complete(t, step_count, *st_n);

        while (t + time_epsilon < t_max) {
            double dt = dt_initial > 0.0 ? dt_initial : timer->compute_dt(*st_n, t);
            if (t + dt > t_max) dt = t_max - t;
            if (dt <= time_epsilon) break;

            st_n = step(t, dt, *st_n, step_count + 1);

            t += dt;
            step_count++;
        }

        for (auto& obs : observers) obs->on_simulation_end();
        return st_n;
    }
};

} // namespace sim

namespace num {
using namespace sim;

/**
 * @brief Compatibility layer for manual assembly pattern.
 * @deprecated Use sim::Simulation inheritance instead.
 */
class SimulationEngine : public Simulation {
public:
    SimulationEngine(std::shared_ptr<geo::IGrid> g,
                     std::shared_ptr<IModel> m,
                     std::shared_ptr<IDiscretizer> d,
                     std::shared_ptr<num::ITimeIntegrator> t,
                     std::shared_ptr<num::ILinearizer> l,
                     std::shared_ptr<num::ISolver> s,
                     std::shared_ptr<utl::IParallelManager> p,
                     std::vector<std::shared_ptr<ISourceSink>> src = {}) {
        this->grd = g;
        this->mdl = m;
        this->discretizer = d;
        this->timer = t;
        this->linearizer = l;
        this->solver = s;
        this->parallel = p;
        this->sources = std::move(src);
        if (this->linearizer) this->linearizer->set_sources(this->sources);
    }

    void build(const utl::ConfigReader& config) override {
        // Already built via constructor
    }
};

} // namespace num
