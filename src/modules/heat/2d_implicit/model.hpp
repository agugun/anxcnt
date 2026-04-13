#pragma once
#include "lib/modules.hpp"
#include "state.hpp"
#include "lib/operators.hpp"

namespace mod {
using namespace top;
namespace heat {

/**
 * @brief 2D Implicit Heat Conduction Model.
 * 
 * Governing Equation:
 *   dT/dt = alpha * (d^2T/dx^2 + d^2T/dy^2)
 * 
 * Physical Logic:
 *   Models planar heat transfer with Dirichlet boundary conditions on all four 
 *   sides (Top, Bottom, Left, Right).
 * 
 * Discretization:
 *   - Spatial: 5-point stencil for the 2D Laplacian.
 *   - Temporal: Backward Euler (Implicit) method.
 * 
 * Boundary Handling:
 *   The boundary cells in the residual are fixed to (T - T_bc), which forces 
 *   the solver to keep them at the prescribed constant temperatures.
 * 
 * Jacobian:
 *   The Jacobian is applied as a linear operator: Jv = (I - dt * alpha * L)v.
 *   Boundary entries in the Jacobian application are simplified to identity 
 *   to maintain Dirichlet constraints during linear solves.
 */
class Heat2DModel : public IModel {
private:
    double alpha;
    double T_top, T_bottom, T_left, T_right;

public:
    Heat2DModel(double alpha, double top, double bottom, double left, double right)
        : alpha(alpha), T_top(top), T_bottom(bottom), T_left(left), T_right(right) {}

    Vector evaluate_rhs(const IState& state) const override {
        const auto& h_state = dynamic_cast<const Heat2DImplicitState&>(state);
        Vector rhs = mop::laplace_2d(h_state.temperatures, (int)h_state.spatial.nx, (int)h_state.spatial.ny, h_state.spatial.dx, h_state.spatial.dy);
        for (auto& v : rhs) v *= alpha;
        return rhs;
    }

    Vector apply_jacobian(const IState& state, const Vector& v, double dt) const override {
        const auto& h_state = dynamic_cast<const Heat2DImplicitState&>(state);
        Vector Jv = mop::laplace_2d(v, (int)h_state.spatial.nx, (int)h_state.spatial.ny, h_state.spatial.dx, h_state.spatial.dy);
        
        // Implicit Euler: (I - dt * J) v
        Vector res(v.size());
        for (size_t i = 0; i < v.size(); ++i) {
            res[i] = v[i] - dt * alpha * Jv[i];
        }

        // Apply Dirichlet BCs for the linear system delta (v is the delta or test vector)
        // For linear systems solver on delta: J_bc * delta = 0
        int nx = (int)h_state.spatial.nx;
        int ny = (int)h_state.spatial.ny;
        for (int i = 0; i < nx; ++i) {
            res[h_state.spatial.idx(i, 0)] = v[h_state.spatial.idx(i, 0)];
            res[h_state.spatial.idx(i, ny - 1)] = v[h_state.spatial.idx(i, ny - 1)];
        }
        for (int j = 0; j < ny; ++j) {
            res[h_state.spatial.idx(0, j)] = v[h_state.spatial.idx(0, j)];
            res[h_state.spatial.idx(nx - 1, j)] = v[h_state.spatial.idx(nx - 1, j)];
        }

        return res;
    }

    Vector build_residual(const IState& s_new, const IState& s_old, double dt) const override {
        const auto& h_new = dynamic_cast<const Heat2DImplicitState&>(s_new);
        const auto& h_old = dynamic_cast<const Heat2DImplicitState&>(s_old);
        
        Vector lap = mop::laplace_2d(h_new.temperatures, (int)h_new.spatial.nx, (int)h_new.spatial.ny, h_new.spatial.dx, h_new.spatial.dy);
        Vector r(h_new.temperatures.size());
        for (size_t i = 0; i < r.size(); ++i) {
            r[i] = h_new.temperatures[i] - h_old.temperatures[i] - dt * alpha * lap[i];
        }

        // Apply Dirichlet BCs to residual (T - T_bc = 0)
        int nx = (int)h_new.spatial.nx;
        int ny = (int)h_new.spatial.ny;
        for (int i = 0; i < nx; ++i) {
            r[h_new.spatial.idx(i, 0)] = h_new.temperatures[h_new.spatial.idx(i, 0)] - T_bottom;
            r[h_new.spatial.idx(i, ny - 1)] = h_new.temperatures[h_new.spatial.idx(i, ny - 1)] - T_top;
        }
        for (int j = 0; j < ny; ++j) {
            r[h_new.spatial.idx(0, j)] = h_new.temperatures[h_new.spatial.idx(0, j)] - T_left;
            r[h_new.spatial.idx(nx - 1, j)] = h_new.temperatures[h_new.spatial.idx(nx - 1, j)] - T_right;
        }

        return r;
    }

    Matrix build_jacobian(const IState&, double) const override { return {}; }
};

} // namespace heat
} // namespace mod
