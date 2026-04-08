#pragma once
#include "../../lib/model/state/field.hpp"
#include "../../lib/numerical_methods/operator/laplacian.hpp"
#include "../../lib/model/i_model.hpp"
#include <stdexcept>

namespace model {

class Heat1DModel : public model::IModel {
public:
    Heat1DModel(double alpha, double left_bc, double right_bc)
        : alpha(alpha), left_bc(left_bc), right_bc(right_bc) {}

    void build_system(double dt, const model::state::IState& istate, model::SystemMatrix& sys) override {
        // Downcast to Field to access dx
        const auto* state_ptr = dynamic_cast<const state::Field*>(&istate);
        if (!state_ptr) {
            throw std::runtime_error("Heat1DModel requires a Field state.");
        }
        const auto& state = *state_ptr;

        int m = state.nx - 2;
        sys.resize(m);

        // 1. Initialize for the time derivative (phi^{n+1}/dt)
        std::fill(sys.a.begin(), sys.a.end(), 0.0);
        std::fill(sys.b.begin(), sys.b.end(), 1.0 / dt);
        std::fill(sys.c.begin(), sys.c.end(), 0.0);
        
        // 2. Set RHS (phi^n / dt)
        for (int i = 0; i < m; ++i) {
            sys.d[i] = state.values[i + 1] / dt;
        }

        // 3. Accumulate Laplacian contribution (-alpha * Laplacian)
        laplace.assembleWithDx(sys.a, sys.b, sys.c, sys.d, state.dx, -alpha);

        // 4. Handle boundary conditions
        sys.d[0] += (alpha / (state.dx * state.dx)) * left_bc;
        sys.d[m - 1] += (alpha / (state.dx * state.dx)) * right_bc;
    }

    void set_boundary_conditions(const std::vector<double>& bcs) override {
        if (bcs.size() >= 2) {
            left_bc = bcs[0];
            right_bc = bcs[1];
        }
    }

    // Keep applying boundaries to state (if needed elsewhere) or adapt logic
    void applyBoundaryConditions(state::Field& state) {
        state.values[0] = left_bc;
        state.values[state.nx - 1] = right_bc;
    }

private:
    double alpha;
    double left_bc;
    double right_bc;
    numerical_methods::operator_unit::Laplacian laplace;
};

} // namespace model
