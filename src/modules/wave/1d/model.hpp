#pragma once
#include "lib/modules.hpp"
#include "state.hpp"
#include "lib/operators.hpp"

namespace mod {
using namespace top;
namespace wave {

/**
 * @brief 1D Wave Propagation Model.
 * 
 * Governing Equation (2nd Order):
 *   d^2u/dt^2 = c^2 * d^2u/dx^2
 * 
 * Mixed First-Order System:
 *   1. du/dt = v
 *   2. dv/dt = c^2 * d^2u/dx^2
 * 
 * Physical Logic:
 *   Models the propagation of waves (e.g., sound or vibration) with speed 'c'. 
 *   The state tracks both displacement (u) and velocity (v).
 * 
 * Discretization:
 *   - Spatial: 2nd Order Central Difference for the displacement Laplacian.
 *   - Temporal: The model provides the RHS for explicit integration. When 
 *     used with a proper integrator (like Leapfrog), it preserves energy.
 */
class Wave1DModel : public IModel {
private:
    double wave_speed;

public:
    Wave1DModel(double c) : wave_speed(c) {}

    Vector evaluate_rhs(const IState& state) const override {
        const auto& w_state = dynamic_cast<const Wave1DState&>(state);
        size_t n = w_state.u.size();
        Vector rhs(2 * n, 0.0);

        // du/dt = v
        for (size_t i = 0; i < n; ++i) {
            rhs[i] = w_state.v[i];
        }

        // dv/dt = c^2 * L(u)
        Vector lap = mop::laplace_1d(w_state.u, w_state.spatial.dx);
        double c2 = wave_speed * wave_speed;
        for (size_t i = 0; i < n; ++i) {
            rhs[n + i] = c2 * lap[i];
        }

        return rhs;
    }

    // Placeholders for implicit methods
    Vector build_residual(const IState&, const IState&, double) const override { return {}; }
    Matrix build_jacobian(const IState&, double) const override { return {}; }
};

} // namespace wave
} // namespace mod
