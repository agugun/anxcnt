/**
 * @file simulation_engine.hpp
 * @brief High-level Orchestrator for Time-Stepped Simulations.
 */
#pragma once
#include "interfaces.hpp"
#include <vector>
#include <memory>

namespace top {

class SimulationEngine {
private:
    std::shared_ptr<IGrid> grid;
    std::shared_ptr<IModel> model;
    std::shared_ptr<IDiscretizer> discretizer;
    std::shared_ptr<ITimeIntegrator> timer;
    std::shared_ptr<ILinearizer> linearizer;
    std::shared_ptr<ISolver> solver;
    std::shared_ptr<IParallelManager> parallel;
    std::vector<std::shared_ptr<ISourceSink>> sources;
    std::vector<std::shared_ptr<IObserver>> observers;

public:
    SimulationEngine(std::shared_ptr<IGrid> g,
                     std::shared_ptr<IModel> m,
                     std::shared_ptr<IDiscretizer> d,
                     std::shared_ptr<ITimeIntegrator> t,
                     std::shared_ptr<ILinearizer> l,
                     std::shared_ptr<ISolver> s,
                     std::shared_ptr<IParallelManager> p,
                     std::vector<std::shared_ptr<ISourceSink>> src = {})
        : grid(g), model(m), discretizer(d), timer(t), linearizer(l), 
          solver(s), parallel(p), sources(std::move(src)) {}

    void add_observer(std::shared_ptr<IObserver> obs) {
        observers.push_back(obs);
    }

    /**
     * @brief Perform a single time step.
     */
    void simulate_step(double dt, IState& state) {
        auto state_next = linearizer->solve_timestep(
            state, dt, *grid, *model, *discretizer, *timer, *solver, *parallel, sources
        );
        state.apply_update(subtract(state_next->to_vector(), state.to_vector()));
    }

    /**
     * @brief Run the full simulation until t_max.
     */
    std::unique_ptr<IState> simulate(double t_max, double dt_initial, std::unique_ptr<IState> initial_state) {
        auto state_n = std::move(initial_state);
        double t = 0.0;
        int step_count = 0;

        for (auto& obs : observers) obs->on_simulation_start(*grid);
        for (auto& obs : observers) obs->on_step_complete(t, step_count, *state_n);

        while (t < t_max) {
            double dt = timer->get_next_timestep(*state_n, t);
            if (t + dt > t_max) dt = t_max - t;

            auto state_next = linearizer->solve_timestep(
                *state_n, dt, *grid, *model, *discretizer, *timer, *solver, *parallel, sources
            );

            state_n = std::move(state_next);
            t += dt;
            step_count++;

            for (auto& obs : observers) obs->on_step_complete(t, step_count, *state_n);
        }

        for (auto& obs : observers) obs->on_simulation_end();
        return state_n;
    }

private:
    Vector subtract(const Vector& a, const Vector& b) {
        Vector res(a.size());
        for (size_t i = 0; i < a.size(); ++i) res[i] = a[i] - b[i];
        return res;
    }
};

} // namespace top
