#pragma once
#include "lib/modules.hpp"
#include "state.hpp"

namespace mod {
using namespace top;
namespace physics_mba {

/**
 * @brief Material Balance Model (Zero-Dimensional / Tank Model).
 * 
 * Governing Equation:
 *   dP/dt = -q / (V * ct)
 * 
 * Physical Logic:
 *   Predicts average reservoir pressure decline over time based on cumulative 
 *   production. It treats the reservoir as a single homogeneous tank.
 *   - V: Total reservoir volume (barrels).
 *   - ct: Total system compressibility (1/psi).
 *   - q: Constant production rate (STB/day).
 * 
 * Numerical Logic:
 *   This is a first-order ODE. While the RHS is constant in this implementation, 
 *   it can be evolved into more complex forms (e.g., pressure-dependent rate).
 * 
 * Assumptions:
 *   - Homogeneous properties.
 *   - Instantaneous pressure equilibrium throughout the tank.
 *   - Single-phase fluid or effective total compressibility.
 */
class MBModel : public IModel {
private:
    double volume; // Reservoir Volume (barrels)
    double ct;     // Total Compressibility (1/psi)
    double q;      // Production Rate (STB/day)

public:
    MBModel(double V, double ct, double q) : volume(V), ct(ct), q(q) {}

    // dP/dt = -q / (V * ct)
    Vector evaluate_rhs(const IState& state) const override {
        Vector rhs(1);
        rhs[0] = -q / (volume * ct);
        return rhs;
    }

    // Placeholders for implicit methods (not strictly needed for this simple linear ODE)
    Vector build_residual(const IState&, const IState&, double) const override { return {}; }
    Matrix build_jacobian(const IState&, double) const override { return {}; }
};

} // namespace physics_mba
} // namespace mod
