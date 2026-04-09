#pragma once
#include "lib/base.hpp"
#include "state.hpp"
#include "lib/operators.hpp"

namespace numerical_methods {
namespace physics_pressure {

class Pressure1DModel : public IModel {
private:
    double eta; // Diffusivity constant (k / (phi * mu * ct))
    double p_left, p_right;

public:
    Pressure1DModel(double k, double phi, double mu, double ct, double pl, double pr)
        : p_left(pl), p_right(pr) {
        eta = 0.0002637 * k / (phi * mu * ct);
    }

    Vector evaluate_rhs(const IState& state) const override {
        const auto& p_state = dynamic_cast<const Pressure1DState&>(state);
        Vector rhs = operators::laplace_1d(p_state.pressures, p_state.dx);
        for (auto& v : rhs) v *= eta;
        return rhs;
    }

    Vector build_residual(const IState& s_new, const IState& s_old, double dt) const override {
        const auto& p_new = dynamic_cast<const Pressure1DState&>(s_new);
        const auto& p_old = dynamic_cast<const Pressure1DState&>(s_old);
        
        Vector lap = operators::laplace_1d(p_new.pressures, p_new.dx);
        Vector r(p_new.pressures.size());
        for (size_t i = 0; i < r.size(); ++i) {
            r[i] = p_new.pressures[i] - p_old.pressures[i] - dt * eta * lap[i];
        }
        r[0] = p_new.pressures[0] - p_left;
        r[r.size()-1] = p_new.pressures[r.size()-1] - p_right;
        return r;
    }

    Matrix build_jacobian(const IState& state, double dt) const override {
        const auto& p_state = dynamic_cast<const Pressure1DState&>(state);
        size_t n = p_state.pressures.size();
        double dx = p_state.dx;
        double coef = eta * dt / (dx * dx);

        Matrix jac(n, Vector(n, 0.0));
        for (size_t i = 1; i < n - 1; ++i) {
            jac[i][i-1] = -coef;
            jac[i][i] = 1.0 + 2.0 * coef;
            jac[i][i+1] = -coef;
        }
        jac[0][0] = 1.0;
        jac[n-1][n-1] = 1.0;
        return jac;
    }
};

} // namespace physics_pressure
} // namespace numerical_methods
