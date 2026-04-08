#include "i_simulator.hpp"
#include "../model/system_matrix.hpp"
#include "../model/state/field.hpp"
#include <iostream>
#include <stdexcept>

namespace simulation {

/**
 * @brief Simple Transient Simulator that wraps around IModel and ISolver.
 *        Compatible with interactive steps (UI) and batch runs.
 */
class TransientSimulator : public ISimulator {
public:
    void set_solver(std::shared_ptr<numerical_methods::solver::ISolver> s) override {
        solver = s;
    }

    /**
     * @brief Run a full batch time loop.
     */
    void run(model::IModel& model, model::state::IState& state, const std::map<std::string, double>& time_params) override {
        double dt = time_params.at("dt");
        double steps = time_params.at("steps");
        double current_time = 0.0;

        for (int i = 0; i < static_cast<int>(steps); ++i) {
            step(dt, model, state);
            current_time += dt;
            // Optionally: write_output(state, current_time);
        }
    }

    /**
     * @brief Execute one explicit implicit timestep.
     */
    void step(double dt, model::IModel& model, model::state::IState& state) {
        if (!solver) throw std::runtime_error("Solver not set in simulator.");
        
        model::SystemMatrix sys;
        model.build_system(dt, state, sys);
        solver->solve(sys, state);
        
        // Post-process (e.g. boundaries are maintained inside build_system typically, 
        // but we can apply them if IModel supports it).
    }

    void write_output(std::ostream& os, const model::state::IState& state, double current_time) override {
        const auto* field = dynamic_cast<const model::state::Field*>(&state);
        if (!field) return;

        // Tidy data format: Time, Index, X, Value
        const auto& data = field->get_data();
        for (int i = 0; i < field->nx; ++i) {
            os << current_time << "," 
               << i << "," 
               << field->get_coordinate(i) << "," 
               << data[i] << "\n";
        }
    }

private:
    std::shared_ptr<numerical_methods::solver::ISolver> solver;
};

} // namespace simulation
