#pragma once
#include "lib/base.hpp"
#include "state.hpp"

namespace numerical_methods {
namespace physics_mba {

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
} // namespace numerical_methods
