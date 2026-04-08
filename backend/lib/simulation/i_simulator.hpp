#pragma once
#include <map>
#include <string>
#include <memory>
#include <ostream>
#include "../model/i_model.hpp"
#include "../model/state/i_state.hpp"
#include "../numerical_methods/solver/i_solver.hpp"

namespace simulation {

/**
 * @brief The Simulator Abstraction (Time Evolution)
 */
class ISimulator {
public:
    virtual ~ISimulator() = default;

    /**
     * @brief Injects the algebraic execution engine (Linear/Non-linear).
     */
    virtual void set_solver(std::shared_ptr<numerical_methods::solver::ISolver> solver) = 0;

    /**
     * @brief The master clock. Loops through time, asks the Model to build the system,
     * asks the Solver to solve it, and updates the State.
     */
    virtual void run(model::IModel& model, model::state::IState& state, const std::map<std::string, double>& time_params) = 0;

    /**
     * @brief Dumps the state to visualization files (or standard output).
     */
    virtual void write_output(std::ostream& os, const model::state::IState& state, double current_time) = 0;
};

} // namespace simulation
