#pragma once
#include "lib/modules.hpp"
#include "state.hpp"
#include "lib/operators.hpp"
#include "modules/reservoir/well.hpp"
#include <cmath>
#include <iostream>
#include <vector>

namespace mod {
using namespace top;
namespace reservoir {

/**
 * @brief 2D Dual-Phase (Oil-Water) Reservoir Model.
 * 
 * Governing Equations:
 *   1. Water: div( (k*krw)/muw * grad(P) ) = d(phi*Sw)/dt
 *   2. Oil:   div( (k*kro)/muo * grad(P) ) = d(phi*So)/dt
 * 
 * Physical Logic:
 *   Models immiscible displacement of oil by water. Both phases are considered 
 *   slightly compressible. Relies on Corey's model for relative permeability.
 * 
 * Discretization:
 *   - Spatial: Finite Volume Method on 5-point Cartesian stencil.
 *   - Upwinding: Single-point upstream weighting for phase mobilities.
 * 
 * Numerical Scheme:
 *   - Fully Implicit (FIM): Both P and Sw are solved simultaneously in the residual.
 *   - Supports an "Implicit Pressure, Explicit Saturation" (IMPES) derivative 
 *     logic through the Jacobian operator if needed.
 * 
 * Stability:
 *   The FIM approach ensures stability even at high flow velocities where 
 *   explicit saturation updates would hit the CFL limit.
 */
class ReservoirDual2DModel : public IModel {
public:
    double k_abs;      // absolute permeability [mD]
    double phi;        // porosity [fraction]
    double mu_w, mu_o; // viscosities [cP]
    double h;          // thickness [ft]
    
    std::vector<std::shared_ptr<ISourceSink>> sources;

    double sw_res;     // residual water saturation
    double so_res;     // residual oil saturation

    ReservoirDual2DModel(double k, double phi_val, double mw, double mo, double h_val, 
                         const std::vector<std::shared_ptr<ISourceSink>>& sources_val)
        : k_abs(k), phi(phi_val), mu_w(mw), mu_o(mo), h(h_val), sources(sources_val), 
          sw_res(0.2), so_res(0.2) {}

    // Relative Permeability (Corey Model)
    void get_rel_perm(double sw, double& krw, double& kro) const {
        double swe = (sw - sw_res) / (1.0 - sw_res - so_res);
        swe = std::max(0.0, std::min(1.0, swe));
        krw = swe * swe;         // nw=2
        kro = (1.0 - swe) * (1.0 - swe); // no=2
    }

    // Fully Implicit Residual for coupled [P, Sw] system
    Vector build_residual(const IState& s_new, const IState& s_old, double dt) const override {
        const auto& state = dynamic_cast<const ReservoirDualPhase2DState&>(s_new);
        const auto& state_old = dynamic_cast<const ReservoirDualPhase2DState&>(s_old);
        int nx = state.spatial.nx, ny = state.spatial.ny;
        double dx = state.spatial.dx, dy = state.spatial.dy;
        int n = (int)state.pressures.size();
        
        Vector res(2 * n, 0.0);
        double unit_conv = 0.001127;
        double pore_vol = phi * (dx * dy * h) / 5.615; // STB

        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                int c = state.idx(i, j);
                int r_w = 2 * c;     // index for Water residual
                int r_o = 2 * c + 1; // index for Oil residual

                // 1. Accumulation Terms
                double comp = 1e-4; // Small compressibility to anchor pressure
                res[r_w] = (pore_vol / dt) * (state.water_saturations[c] - state_old.water_saturations[c]) 
                           + comp * (state.pressures[c] - state_old.pressures[c]);
                res[r_o] = (pore_vol / dt) * ((1.0 - state.water_saturations[c]) - (1.0 - state_old.water_saturations[c]));

                // 2. Flux Terms
                auto add_fluxes = [&](int ni, int nj, double dist) {
                    int n_idx = state.idx(ni, nj);
                    double p_c = state.pressures[c];
                    double p_n = state.pressures[n_idx];
                    
                    // Upwinding
                    int up = (p_c > p_n) ? c : n_idx;
                    double krw, kro;
                    get_rel_perm(state.water_saturations[up], krw, kro);
                    
                    double Tw = unit_conv * k_abs * (krw / mu_w) * (h * (dist == dx ? dy : dx)) / dist;
                    double To = unit_conv * k_abs * (kro / mu_o) * (h * (dist == dx ? dy : dx)) / dist;
                    
                    res[r_w] -= Tw * (p_n - p_c);
                    res[r_o] -= To * (p_n - p_c);
                };

                if (i > 0) add_fluxes(i - 1, j, dx);
                if (i < nx - 1) add_fluxes(i + 1, j, dx);
                if (j > 0) add_fluxes(i, j - 1, dy);
                if (j < ny - 1) add_fluxes(i, j + 1, dy);
            }
        }

        // 3. Source Terms (Abstracted)
        for (auto& s : sources) {
            s->apply(res, nullptr, state, dt);
        }

        return res;
    }

    Matrix build_jacobian(const IState& state_raw, double dt) const override {
        const auto& r_state = dynamic_cast<const ReservoirDualPhase2DState&>(state_raw);
        int nx = r_state.spatial.nx, ny = r_state.spatial.ny;
        double dx = r_state.spatial.dx, dy = r_state.spatial.dy;
        int n = (int)r_state.pressures.size();
        
        Matrix jac(2 * n, Vector(2 * n, 0.0));
        double unit_conv = 0.001127;
        double pore_vol = phi * (dx * dy * h) / 5.615;
        double den_sw = 1.0 - sw_res - so_res;

        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                int c = r_state.idx(i, j);
                int r_w = 2 * c;
                int r_o = 2 * c + 1;

                // 1. Accumulation Jacobian
                double comp = 1e-4;
                jac[r_w][2 * c] = comp;            // dRw/dP
                jac[r_w][r_w + 1] = pore_vol / dt;  // dRw/dSw
                jac[r_o][r_w + 1] = -pore_vol / dt; // dRo/dSw

                // 2. Flux Jacobian
                auto add_flux_jac = [&](int ni, int nj, double dist) {
                    int n_idx = r_state.idx(ni, nj);
                    int n_p = 2 * n_idx;
                    int n_s = 2 * n_idx + 1;
                    
                    double p_c = r_state.pressures[c];
                    double p_n = r_state.pressures[n_idx];
                    double dp = p_n - p_c;
                    
                    // Upwinding
                    int up = (p_c > p_n) ? c : n_idx;
                    double sw_up = r_state.water_saturations[up];
                    double krw, kro;
                    get_rel_perm(sw_up, krw, kro);
                    
                    // RelPerm Derivatives (Continuity at Sw limits)
                    double dkrw_dsw = 0.0;
                    double dkro_dsw = 0.0;
                    if (sw_up > sw_res && sw_up < (1.0 - so_res)) {
                        double swe = (sw_up - sw_res) / den_sw;
                        dkrw_dsw = 2.0 * swe / den_sw;
                        dkro_dsw = -2.0 * (1.0 - swe) / den_sw;
                    }

                    double Tw = unit_conv * k_abs * (krw / mu_w) * (h * (dist == dx ? dy : dx)) / dist;
                    double To = unit_conv * k_abs * (kro / mu_o) * (h * (dist == dx ? dy : dx)) / dist;
                    double dTw_dsw = unit_conv * k_abs * (dkrw_dsw / mu_w) * (h * (dist == dx ? dy : dx)) / dist;
                    double dTo_dsw = unit_conv * k_abs * (dkro_dsw / mu_o) * (h * (dist == dx ? dy : dx)) / dist;

                    // dRw/dPi, dRw/dPj, dRo/dPi, dRo/dPj
                    jac[r_w][2 * c] += Tw;
                    jac[r_w][n_p]   -= Tw;
                    jac[r_o][2 * c] += To;
                    jac[r_o][n_p]   -= To;

                    // dR/dSw_upwind
                    if (up == c) {
                        jac[r_w][2 * c + 1] -= dTw_dsw * dp;
                        jac[r_o][2 * c + 1] -= dTo_dsw * dp;
                    } else {
                        jac[r_w][n_s]       -= dTw_dsw * dp;
                        jac[r_o][n_s]       -= dTo_dsw * dp;
                    }
                };

                if (i > 0) add_flux_jac(i - 1, j, dx);
                if (i < nx - 1) add_flux_jac(i + 1, j, dx);
                if (j > 0) add_flux_jac(i, j - 1, dy);
                if (j < ny - 1) add_flux_jac(i, j + 1, dy);
            }
        }

        // 3. Source Jacobian (Abstracted)
        Vector dummy_res(2 * n, 0.0);
        for (auto& s : sources) {
            s->apply(dummy_res, &jac, state_raw, dt);
        }

        return jac;
    }

    Vector apply_jacobian(const IState& state_raw, const Vector& v, double dt) const override {
        const auto& state = dynamic_cast<const ReservoirDualPhase2DState&>(state_raw);
        int nx = state.spatial.nx, ny = state.spatial.ny;
        double dx = state.spatial.dx, dy = state.spatial.dy;
        Vector res(v.size(), 0.0);
        double unit_conv = 0.001127;

        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                int cur = state.idx(i, j);
                double sw = state.water_saturations[cur];
                double krw, kro;
                get_rel_perm(sw, krw, kro);
                double lambda_t = (krw / mu_w) + (kro / mu_o);

                double diag = 0.0;
                auto apply_neighbor = [&](int ni, int nj, double dist) {
                    int neighbor = state.idx(ni, nj);
                    double sw_n = state.water_saturations[neighbor];
                    double krw_n, kro_n;
                    get_rel_perm(sw_n, krw_n, kro_n);
                    double lambda_t_n = (krw_n / mu_w) + (kro_n / mu_o);
                    
                    double lambda_avg = 0.5 * (lambda_t + lambda_t_n);
                    double T = unit_conv * k_abs * lambda_avg * (h * (dist == dx ? dy : dx)) / dist;
                    
                    res[cur] += T * v[neighbor];
                    diag -= T;
                };

                if (i > 0) apply_neighbor(i - 1, j, dx);
                if (i < nx - 1) apply_neighbor(i + 1, j, dx);
                if (j > 0) apply_neighbor(i, j - 1, dy);
                if (j < ny - 1) apply_neighbor(i, j + 1, dy);

                res[cur] += diag * v[cur];
            }
        }
        return res;
    }

    // Evaluate RHS (not used in IMPES exactly like this, but required by interface)
    Vector evaluate_rhs(const IState& state) const override { return {}; }
};

} // namespace reservoir
} // namespace mod
