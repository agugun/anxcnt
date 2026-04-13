#include "lib/modules.hpp"
#include "state.hpp"
#include "lib/operators.hpp"

namespace mod {
using namespace top;
namespace heat {

/**
 * @brief 1D Explicit Heat Conduction Model.
 * 
 * Governing Equation:
 *   dT/dt = alpha * d^2T/dx^2
 * 
 * Physical Logic:
 *   The model calculates the rate of temperature change based on Fourier's Law
 *   of heat conduction. alpha represents the thermal diffusivity (k / (rho * Cp)).
 * 
 * Discretization:
 *   - Spatial: 2nd Order Central Difference scheme (mop::laplace_1d).
 *   - Temporal: The RHS (Right Hand Side) is evaluated at the current time step (n)
 *     and is meant to be integrated using Explicit schemes like Forward Euler.
 * 
 * Numerical Stability:
 *   Requires CFL condition: dt <= dx^2 / (2 * alpha).
 */
class HeatModel : public IModel {
private:
    double alpha;
    Spatial1D spatial;

public:
    HeatModel(double alpha, Spatial1D spatial) : alpha(alpha), spatial(spatial) {}

    /**
     * @brief Computes the Right-Hand-Side of the heat equation.
     * @return alpha * (T[i+1] - 2T[i] + T[i-1]) / dx^2
     */
    Vector evaluate_rhs(const IState& state) const override {
        const auto& heat_state = dynamic_cast<const Heat1DExplicitState&>(state);
        Vector rhs = mop::laplace_1d(heat_state.temperatures, spatial.dx);
        
        // Multiply by alpha: scale the Laplacian by diffusivity
        for (auto& val : rhs) val *= alpha;
        
        return rhs;
    }


    // Dummy implementations for implicit methods
    Vector build_residual(const IState&, const IState&, double) const override { return {}; }
    Matrix build_jacobian(const IState&, double) const override { return {}; }
};

} // namespace heat
} // namespace mod