#pragma once
#include "lib/modules.hpp"
#include "state.hpp"
#include "lib/operators.hpp"

namespace mod {
using namespace top;
namespace physics_heat {

class Heat3DModel : public IModel {
private:
    double alpha;
    double T_front, T_back, T_top, T_bottom, T_left, T_right;

public:
    Heat3DModel(double alpha, double front, double back, double top, double bottom, double left, double right)
        : alpha(alpha), T_front(front), T_back(back), T_top(top), T_bottom(bottom), T_left(left), T_right(right) {}

    Vector evaluate_rhs(const IState& state) const override {
        const auto& h_state = dynamic_cast<const Heat3DState&>(state);
        Vector rhs = mop::laplace_3d(h_state.temperatures, h_state.nx, h_state.ny, h_state.nz, h_state.dx, h_state.dy, h_state.dz);
        for (auto& v : rhs) v *= alpha;
        return rhs;
    }

    Vector apply_jacobian(const IState& state, const Vector& v, double dt) const override {
        const auto& h_state = dynamic_cast<const Heat3DState&>(state);
        Vector Jv = mop::laplace_3d(v, h_state.nx, h_state.ny, h_state.nz, h_state.dx, h_state.dy, h_state.dz);
        
        // Implicit Euler: (I - dt * J) v
        Vector res(v.size());
        for (size_t i = 0; i < v.size(); ++i) {
            res[i] = v[i] - dt * alpha * Jv[i];
        }

        // Apply Dirichlet BCs for the linear system (boundary deltas should be zero or specifically handled)
        int nx = h_state.nx;
        int ny = h_state.ny;
        int nz = h_state.nz;

        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                res[h_state.idx(0, j, k)] = v[h_state.idx(0, j, k)];
                res[h_state.idx(nx - 1, j, k)] = v[h_state.idx(nx - 1, j, k)];
            }
        }
        for (int k = 0; k < nz; ++k) {
            for (int i = 0; i < nx; ++i) {
                res[h_state.idx(i, 0, k)] = v[h_state.idx(i, 0, k)];
                res[h_state.idx(i, ny - 1, k)] = v[h_state.idx(i, ny - 1, k)];
            }
        }
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                res[h_state.idx(i, j, 0)] = v[h_state.idx(i, j, 0)];
                res[h_state.idx(i, j, nz - 1)] = v[h_state.idx(i, j, nz - 1)];
            }
        }

        return res;
    }

    Vector build_residual(const IState& s_new, const IState& s_old, double dt) const override {
        const auto& h_new = dynamic_cast<const Heat3DState&>(s_new);
        const auto& h_old = dynamic_cast<const Heat3DState&>(s_old);
        
        Vector lap = mop::laplace_3d(h_new.temperatures, h_new.nx, h_new.ny, h_new.nz, h_new.dx, h_new.dy, h_new.dz);
        Vector r(h_new.temperatures.size());
        for (size_t i = 0; i < r.size(); ++i) {
            r[i] = h_new.temperatures[i] - h_old.temperatures[i] - dt * alpha * lap[i];
        }

        // Apply Dirichlet BCs to residual (T - T_bc = 0)
        int nx = h_new.nx;
        int ny = h_new.ny;
        int nz = h_new.nz;

        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                r[h_new.idx(0, j, k)] = h_new.temperatures[h_new.idx(0, j, k)] - T_left;
                r[h_new.idx(nx - 1, j, k)] = h_new.temperatures[h_new.idx(nx - 1, j, k)] - T_right;
            }
        }
        for (int k = 0; k < nz; ++k) {
            for (int i = 0; i < nx; ++i) {
                r[h_new.idx(i, 0, k)] = h_new.temperatures[h_new.idx(i, 0, k)] - T_bottom;
                r[h_new.idx(i, ny - 1, k)] = h_new.temperatures[h_new.idx(i, ny - 1, k)] - T_top;
            }
        }
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                r[h_new.idx(i, j, 0)] = h_new.temperatures[h_new.idx(i, j, 0)] - T_front;
                r[h_new.idx(i, j, nz - 1)] = h_new.temperatures[h_new.idx(i, j, nz - 1)] - T_back;
            }
        }

        return r;
    }

    Matrix build_jacobian(const IState&, double) const override { return {}; }
};

} // namespace physics_heat
} // namespace mod
