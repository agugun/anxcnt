#pragma once
#include "../../model/system_matrix.hpp"
#include "../../model/state/i_state.hpp"

namespace numerical_methods::solver {

/**
 * @brief Represents the algebraic execution engine (Linear/Non-linear).
 */
class ISolver {
public:
    virtual ~ISolver() = default;

    /**
     * @brief Solves the assembled system and updates the state.
     */
    virtual void solve(const model::SystemMatrix& sys, model::state::IState& state) = 0;
};

} // namespace numerical_methods::solver
