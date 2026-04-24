#pragma once
#include "lib/interfaces.hpp"
#include "lib/integrators.hpp"
#include "lib/simulation_engine.hpp"
#include "model.hpp"
#include "state.hpp"

namespace mod::oscillator {

/**
 * @brief Discretizer for the Harmonic Oscillator.
 * Translates the analytical model rules into the Residual and Jacobian matrices.
 * This represents Phase 3: The C++ Numerics (Flexibility & Semantics).
 */
class OscillatorDiscretizer : public top::IDiscretizer {
public:
    void assemble_jacobian(const top::IGrid& grid, const top::IModel& base_model, const top::IState& base_state, top::SparseMatrix& J) const override {
        const auto& model = dynamic_cast<const OscillatorModel&>(base_model);
        
        // J = dR/du
        // Eq. 1: dx/dt - v = 0      --> dR0/dx = 0, dR0/dv = -1
        // Eq. 2: m*dv/dt + cx + kv = 0 --> dR1/dx = k, dR1/dv = c
        // (Note: The Accumulation terms dR_acc/du are automatically added by the Integrator)

        J.triplets.push_back({0, 1, -1.0});         // dR0 / dv
        J.triplets.push_back({1, 0, model.k});      // dR1 / dx
        J.triplets.push_back({1, 1, model.c});      // dR1 / dv
    }

    void assemble_residual(const top::IGrid& grid, const top::IModel& base_model, const top::IState& base_state, top::Vector& R) const override {
        const auto& model = dynamic_cast<const OscillatorModel&>(base_model);
        const auto& state = dynamic_cast<const OscillatorState&>(base_state);

        // Vectorized Math Semantics: The math reads exactly like the equations.
        // Eq. 1 [Notebook] Physics Residual for Position
        double r_pos = -state.v; 
        
        // Eq. 2 [Notebook] Physics Residual for Velocity 
        double r_vel = model.k * state.x + model.c * state.v;

        R[0] += r_pos;
        R[1] += r_vel;
    }

    void apply_boundary_conditions(const top::IGrid& grid, const top::IModel& base_model, const top::IState& base_state, top::SparseMatrix& J, top::Vector& R) const override {
        // No boundaries for a 0D ODE problem.
    }
};

/**
 * @brief Dummy grid for the 2-variable ODE system.
 */
class OscillatorGrid : public top::IGrid {
public:
    size_t get_total_cells() const override { return 2; }
};

} // namespace mod::oscillator
