#pragma once
#include "../../lib/model/state/field.hpp"
#include "../../lib/numerical_methods/operator/laplacian.hpp"
#include "../../lib/model/i_model.hpp"
#include <stdexcept>

namespace model {

class Pressure1DModel : public model::IModel {
public:
    Pressure1DModel(double k, double phi, double mu, double ct, double left_bc, double right_bc)
        : k(k), phi(phi), mu(mu), ct(ct), left_bc(left_bc), right_bc(right_bc), p_prev(nullptr) {}

    void set_previous_state(const model::state::IState& prev_state) override {
        p_prev = dynamic_cast<const state::Field*>(&prev_state);
        if (!p_prev) {
            throw std::runtime_error("Pressure1DModel requires a Field state for p_prev.");
        }
    }

    void build_system(double dt, const model::state::IState& istate, model::SystemMatrix& sys) override {
        const auto* state_ptr = dynamic_cast<const state::Field*>(&istate);
        if (!state_ptr) {
            throw std::runtime_error("Pressure1DModel requires a Field state.");
        }
        if (!p_prev) {
            throw std::runtime_error("Pressure1DModel requires p_prev to be set before building the system.");
        }

        const auto& p_curr = *state_ptr;
        int m = p_curr.nx - 2;
        sys.resize(m);

        double eta = 0.0002637 * k / (phi * mu * ct);

        // 1. Initialize Jacobian with Identity derivative (1/dt)
        std::fill(sys.a.begin(), sys.a.end(), 0.0);
        std::fill(sys.b.begin(), sys.b.end(), 1.0 / dt);
        std::fill(sys.c.begin(), sys.c.end(), 0.0);
        
        // 2. Initialize Residual with (p_curr - p_prev) / dt
        // In our SystemMatrix, 'd' stores the RHS / Residual vector
        double inv_dt = 1.0 / dt;
        for (int i = 0; i < m; ++i) {
            sys.d[i] = (p_curr.values[i + 1] - p_prev->values[i + 1]) * inv_dt;
        }

        // 3. Accumulate Laplacian contribution to Jacobian (matrix)
        laplace.assembleWithDx(sys.a, sys.b, sys.c, sys.d, p_curr.dx, -eta);

        // 4. Accumulate Laplacian contribution to Residual (vector)
        // R = (p - p_prev)/dt - eta * L(p)
        double dx2 = p_curr.dx * p_curr.dx;
        for (int i = 0; i < m; ++i) {
            int idx = i + 1;
            double p_l = p_curr.values[idx - 1];
            double p_m = p_curr.values[idx];
            double p_r = p_curr.values[idx + 1];
            sys.d[i] -= eta * (p_l - 2.0 * p_m + p_r) / dx2;
        }
    }

    void set_boundary_conditions(const std::vector<double>& bcs) override {
        if (bcs.size() >= 2) {
            left_bc = bcs[0];
            right_bc = bcs[1];
        }
    }

    // Keep original method for where it's needed explicitly
    void applyBoundaryConditions(state::Field &state) {
        state.values[0] = left_bc;
        state.values[state.nx - 1] = right_bc;
    }

private:
    double k, phi, mu, ct;
    double left_bc, right_bc;
    numerical_methods::operator_unit::Laplacian laplace;
    const state::Field* p_prev;
};

} // namespace model
