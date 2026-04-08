#pragma once
#include <vector>
#include "state/i_state.hpp"
#include "system_matrix.hpp"

namespace model {

/**
 * @brief The Model/Space Abstraction (Physics)
 */
class IModel {
public:
    virtual ~IModel() = default;

    /**
     * @brief Defines how the edges of the spatial domain behave.
     */
    virtual void set_boundary_conditions(const std::vector<double>& bcs) = 0;

    /**
     * @brief Applies the governing equations and spatial operators to the State.
     * Returns the assembled matrices/vectors for the solver.
     */
    virtual void build_system(double dt, const state::IState& state, SystemMatrix& sys) = 0;

    /**
     * @brief Optional: Set previous state for transient/nonlinear systems.
     */
    virtual void set_previous_state(const state::IState& prev_state) {}
};

} // namespace model
