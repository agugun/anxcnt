#pragma once
#include "lib/modules.hpp"
#include "state.hpp"
#include "lib/operators.hpp"
#include "modules/reservoir/well.hpp"
#include <cmath>
#include <vector>

namespace mod {
using namespace top;
namespace reservoir {

/**
 * @brief Classical 2D Pressure Diffusivity Reservoir Model.
 * 
 * Governing Equation:
 *   d^2P/dx^2 + d^2P/dy^2 = (1 / eta) * dP/dt
 * 
 * Physical Logic:
 *   Models single-phase fluid flow in a 2D horizontal reservoir layer. 
 *   Transmissibility is assumed to be isotropic unless dx != dy.
 * 
 * Discretization:
 *   - Spatial: Finite Volume Method on a 5-point Cartesian stencil.
 *   - Temporal: Backward Euler (Implicit) method.
 * 
 * Boundary Conditions:
 *   Implements No-flow (Neumann) boundaries on all four sides of the 
 *   rectangular grid by mirroring pressure values at the boundaries.
 * 
 * Jacobian:
 *   Implemented as both a sparse matrix constructor and a linear operator 
 *   (apply_jacobian) for iterative solvers.
 */
class Reservoir2DModel : public IModel {
private:
    double k;      // permeability [mD]
    double phi;    // porosity [fraction]
    double mu;     // viscosity [cP]
    double ct;     // total compressibility [psi^-1]
    double B;      // formation volume factor [rb/stb]
    double h;      // thickness [ft]
    
    std::vector<std::shared_ptr<ISourceSink>> sources;

    double diffusivity; 

public:
    Reservoir2DModel(double k_val, double phi_val, double mu_val, double ct_val, 
                     double B_val, double h_val, const std::vector<std::shared_ptr<ISourceSink>>& sources_val)
        : k(k_val), phi(phi_val), mu(mu_val), ct(ct_val), 
          B(B_val), h(h_val), sources(sources_val) {
        
        diffusivity = 0.0002637 * k / (phi * mu * ct);
    }

    Vector evaluate_rhs(const IState& state) const override {
        const auto& r_state = dynamic_cast<const Reservoir2DState&>(state);
        Vector rhs = mop::laplace_2d(r_state.pressures, r_state.spatial.dx, r_state.spatial.dy, r_state.spatial.nx, r_state.spatial.ny);
        for (auto& v : rhs) v *= diffusivity;

        for (auto& s : sources) {
            s->apply(rhs, nullptr, r_state, 0.0);
        }

        return rhs;
    }

    Vector build_residual(const IState& s_new, const IState& s_old, double dt) const override {
        const auto& r_new = dynamic_cast<const Reservoir2DState&>(s_new);
        const auto& r_old = dynamic_cast<const Reservoir2DState&>(s_old);
        int nx = r_new.spatial.nx, ny = r_new.spatial.ny;
        double dx = r_new.spatial.dx, dy = r_new.spatial.dy;
        
        Vector res(r_new.pressures.size());
        double tx = diffusivity * dt / (dx * dx);
        double ty = diffusivity * dt / (dy * dy);

        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                int cur = r_new.idx(i, j);
                
                // No-flow (Neumann) BCs: mirror value at boundaries
                double p_cur = r_new.pressures[cur];
                double p_l = (i == 0) ? p_cur : r_new.pressures[r_new.idx(i-1, j)];
                double p_r = (i == nx - 1) ? p_cur : r_new.pressures[r_new.idx(i+1, j)];
                double p_b = (j == 0) ? p_cur : r_new.pressures[r_new.idx(i, j-1)];
                double p_t = (j == ny - 1) ? p_cur : r_new.pressures[r_new.idx(i, j+1)];
                
                double lap_p = (p_r - 2.0 * p_cur + p_l) * tx + (p_t - 2.0 * p_cur + p_b) * ty;
                res[cur] = p_cur - r_old.pressures[cur] - lap_p;
            }
        }

        Vector well_contributions(res.size(), 0.0);
        for (auto& s : sources) {
            s->apply(well_contributions, nullptr, r_new, dt);
        }
        for (size_t i = 0; i < res.size(); ++i) {
            res[i] += dt * well_contributions[i];
        }

        return res;
    }

    Matrix build_jacobian(const IState& state, double dt) const override {
        const auto& r_state = dynamic_cast<const Reservoir2DState&>(state);
        int nx = r_state.spatial.nx, ny = r_state.spatial.ny;
        double dx = r_state.spatial.dx, dy = r_state.spatial.dy;
        int n = nx * ny;
        
        double tx = diffusivity * dt / (dx * dx);
        double ty = diffusivity * dt / (dy * dy);

        Matrix jac(n, Vector(n, 0.0));
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                int cur = r_state.idx(i, j);
                
                jac[cur][cur] = 1.0 + 2.0 * tx + 2.0 * ty;
                
                if (i > 0) jac[cur][r_state.idx(i-1, j)] = -tx;
                else jac[cur][cur] -= tx;
                
                if (i < nx - 1) jac[cur][r_state.idx(i+1, j)] = -tx;
                else jac[cur][cur] -= tx;
                
                if (j > 0) jac[cur][r_state.idx(i, j-1)] = -ty;
                else jac[cur][cur] -= ty;
                
                if (j < ny - 1) jac[cur][r_state.idx(i, j+1)] = -ty;
                else jac[cur][cur] -= ty;
            }
        }
        return jac;
    }

    Vector apply_jacobian(const IState& state, const Vector& v, double dt) const override {
        const auto& r_state = dynamic_cast<const Reservoir2DState&>(state);
        int nx = r_state.spatial.nx;
        int ny = r_state.spatial.ny;
        double dx = r_state.spatial.dx;
        double dy = r_state.spatial.dy;
        Vector res(v.size());
        
        double tx = diffusivity * dt / (dx * dx);
        double ty = diffusivity * dt / (dy * dy);

        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                int cur = r_state.idx(i, j);
                
                // (I - dt * eta * L) * v
                double v_cur = v[cur];
                double v_l = (i == 0) ? v_cur : v[r_state.idx(i-1, j)];
                double v_r = (i == nx - 1) ? v_cur : v[r_state.idx(i+1, j)];
                double v_b = (j == 0) ? v_cur : v[r_state.idx(i, j-1)];
                double v_t = (j == ny - 1) ? v_cur : v[r_state.idx(i, j+1)];
                
                double d2v_dx2_dt_eta = (v_r - 2.0 * v_cur + v_l) * tx;
                double d2v_dy2_dt_eta = (v_t - 2.0 * v_cur + v_b) * ty;
                
                res[cur] = v_cur - (d2v_dx2_dt_eta + d2v_dy2_dt_eta);
            }
        }
        return res;
    }
};

} // namespace reservoir
} // namespace mod
