#pragma once
#include "lib/modules.hpp"
#include "state.hpp"
#include "lib/operators.hpp"
#include "lib/spatial.hpp"
 
namespace mod {
using namespace top;
namespace heat {

/**
 * @brief 1D Implicit Heat Conduction Model.
 * 
 * Governing Equation:
 *   dT/dt = alpha * d^2T/dx^2
 * 
 * Physical Logic:
 *   Models heat transfer along a 1D rod. Boundary conditions are Dirichlet
 *   at the left and right ends.
 * 
 * Discretization:
 *   - Spatial: 2nd Order Central Difference scheme for the Laplacian.
 *   - Temporal: Backward Euler (Implicit) method.
 * 
 * Residual Form (F = 0):
 *   F_i = T_i^{n+1} - T_i^n - dt * alpha * (T_{i+1}^{n+1} - 2T_i^{n+1} + T_{i-1}^{n+1}) / dx^2
 * 
 * Jacobian (J = dF/dT):
 *   Diagonal: 1 + 2 * alpha * dt / dx^2
 *   Off-diagonals: -alpha * dt / dx^2
 * 
 * The system is strictly diagonally dominant, ensuring stable convergence.
 */
class Heat1DModel : public IModel {
private:
    double alpha;
    double T_left, T_right;

public:
    Heat1DModel(double alpha, double TL, double TR) 
        : alpha(alpha), T_left(TL), T_right(TR) {}

    Vector evaluate_rhs(const IState& state) const override {
        const auto& h_state = dynamic_cast<const Heat1DImplicitState&>(state);
        Vector rhs = mop::laplace_1d(h_state.temperatures, h_state.spatial.dx);
        for (auto& v : rhs) v *= alpha;
        return rhs;
    }

    Vector apply_jacobian(const IState& state, const Vector& v, double dt) const override {
        const auto& h_state = dynamic_cast<const Heat1DImplicitState&>(state);
        Vector Jv = mop::laplace_1d(v, h_state.spatial.dx);
        // Implicit Euler Jacobian application: (I - dt * J) v
        Vector res(v.size());
        for (size_t i = 1; i < v.size() - 1; ++i) {
            res[i] = v[i] - dt * alpha * Jv[i];
        }
        res[0] = v[0];
        res[v.size()-1] = v[v.size()-1];
        return res;
    }

    Vector build_residual(const IState& s_new, const IState& s_old, double dt) const override {
        const auto& h_new = dynamic_cast<const Heat1DImplicitState&>(s_new);
        const auto& h_old = dynamic_cast<const Heat1DImplicitState&>(s_old);
        
        Vector lap = mop::laplace_1d(h_new.temperatures, h_new.spatial.dx);
        Vector r(h_new.temperatures.size());
        for (size_t i = 1; i < r.size() - 1; ++i) {
            r[i] = h_new.temperatures[i] - h_old.temperatures[i] - dt * alpha * lap[i];
        }
        r[0] = h_new.temperatures[0] - T_left;
        r[r.size()-1] = h_new.temperatures[r.size()-1] - T_right;
        return r;
    }

    Matrix build_jacobian(const IState& state, double dt) const override {
        const auto& h_state = dynamic_cast<const Heat1DImplicitState&>(state);
        size_t n = h_state.spatial.nx;
        double dx = h_state.spatial.dx;
        double coef = alpha * dt / (dx * dx);

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

    void set_bcs(double left, double right) {
        T_left = left;
        T_right = right;
    }
};

} // namespace heat
} // namespace mod
