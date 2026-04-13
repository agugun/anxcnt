#pragma once
#include "lib/modules.hpp"
#include "state.hpp"
#include "lib/operators.hpp"
#include "lib/operators.hpp"
#include <cmath>
#include <vector>

namespace mod {
using namespace top;
namespace reservoir {

class ReservoirWellOilGas2D : public ISourceSink {
public:
    int i, j;
    double q_total;
    bool is_injector;
    
    // Callbacks to access model properties
    std::function<void(double, double&, double&)> get_rel_perm;
    double bo, mu_o, mu_g;
    std::function<double(double)> get_bg;

    ReservoirWellOilGas2D(int i_v, int j_v, double q, bool inj, 
                          std::function<void(double, double&, double&)> rp, 
                          double b_o, double mo, double mg, 
                          std::function<double(double)> g_bg)
        : i(i_v), j(j_v), q_total(q), is_injector(inj), 
          get_rel_perm(rp), bo(b_o), mu_o(mo), mu_g(mg), get_bg(g_bg) {}

    void apply(Vector& residual, Matrix* jacobian, const top::IState& state_raw, double dt,
               std::vector<SparseMatrix::Entry>* sparse_entries = nullptr) override {
        const auto& state = dynamic_cast<const ReservoirOilGas2DState&>(state_raw);
        int c = state.idx(i, j);
        double sw_val = state.gas_saturations[c]; 
        double p_val = state.pressures[c];
        double bg = get_bg(p_val);

        if (is_injector) {
            residual[2 * c + 1] -= q_total; 
            // Jacobian entry for injection: dRes/dSg is 0 (rate is constant)
        } else {
            double krog, krg;
            get_rel_perm(sw_val, krog, krg);
            
            double lam_o = krog / (mu_o * bo);
            double lam_g = krg / (mu_g * bg);
            double lam_t = lam_o + lam_g;
            double fg = lam_g / lam_t;
            
            residual[2 * c] += q_total * (1.0 - fg);
            residual[2 * c + 1] += q_total * fg;

            // Optional Jacobian contributions for wells can be added here
            // but usually they are small compared to flux.
        }
    }
};

/**
 * @brief 2D Oil-Gas Reservoir Model.
 * 
 * Governing Equations:
 *   1. Oil (dead): div( (k*kro)/(muo*Bo) * grad(P) ) = d(phi*So/Bo)/dt
 *   2. Gas:        div( (k*krg)/(mug*Bg(P)) * grad(P) ) = d(phi*Sg/Bg(P))/dt
 * 
 * Physical Logic:
 *   Models displacement of oil by gas or gas production. Assumes "dead oil" 
 *   (constant Bo) and compressible gas following the Real Gas Law (Bg inversely 
 *   proportional to Pressure).
 * Discretization:
 *   - Spatial: Finite Volume on 5-point Cartesian stencil.
 *   - Upwinding: Phase mobilities (mobility/viscosity/formation volume factor) 
 *   - evaluated using single-point upstream weighting.
 * 
 * Numerical Scheme:
 *   - Fully Implicit (FIM) with Newton-Raphson.
 */
class ReservoirOilGas2DModel : public IModel {
public:
    double k_abs;      // absolute permeability [mD]
    double phi;        // porosity [fraction]
    double mu_o, mu_g; // viscosities [cP]
    double h;          // thickness [ft]
    
    // Properties
    double bo;         // Oil formation volume factor (assumed constant for dead oil)
    double p_std;      // Standard pressure [psi]
    double bg_std;     // Gas formation volume factor at standard conditions

    std::vector<std::shared_ptr<ISourceSink>> sources;

    ReservoirOilGas2DModel(double k, double phi_val, double mo, double mg, double h_val, 
                           const std::vector<std::shared_ptr<ISourceSink>>& sources_val)
        : k_abs(k), phi(phi_val), mu_o(mo), mu_g(mg), h(h_val), sources(sources_val), 
          bo(1.2), p_std(14.7), bg_std(1.0) {}

    // Gas Formation Volume Factor: Bg = Bg_std * (P_std / P)
    double get_bg(double p) const {
        return bg_std * (p_std / std::max(1.0, p));
    }

    // Derivative: dBg/dP = -Bg_std * P_std / P^2
    double get_dbg_dp(double p) const {
        return -bg_std * p_std / (std::max(1.0, p) * std::max(1.0, p));
    }

    // Relative Permeability (Gas-Oil Corey)
    void get_rel_perm(double sg, double& krog, double& krg) const {
        double swc = 0.2;
        double sorg = 0.1;
        double sge = sg / (1.0 - swc - sorg);
        sge = std::max(0.0, std::min(1.0, sge));
        krg = sge * sge;
        krog = (1.0 - sge) * (1.0 - sge);
    }

    Vector build_residual(const IState& s_new, const IState& s_old, double dt) const override {
        const auto& state = dynamic_cast<const ReservoirOilGas2DState&>(s_new);
        const auto& state_old = dynamic_cast<const ReservoirOilGas2DState&>(s_old);
        int nx = state.spatial.nx, ny = state.spatial.ny;
        double dx = state.spatial.dx, dy = state.spatial.dy;
        int n = (int)state.pressures.size();
        
        Vector res(2 * n, 0.0);
        double unit_conv = 0.001127;
        double pore_vol_unit = (dx * dy * h) / 5.615;

        #pragma omp parallel for collapse(2)
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                int c = state.idx(i, j);
                double p_c = state.pressures[c];
                double sg_c = state.gas_saturations[c];
                double so_c = 1.0 - sg_c - state.swc;

                double bg_c = get_bg(p_c);

                // 1. Accumulation
                res[2 * c] = (phi * pore_vol_unit / dt) * (so_c / bo - (1.0 - state_old.gas_saturations[c] - state_old.swc) / bo);
                res[2 * c + 1] = (phi * pore_vol_unit / dt) * (sg_c / bg_c - state_old.gas_saturations[c] / get_bg(state_old.pressures[c]));

                // 2. Flux
                auto add_fluxes = [&](int ni, int nj, double dist) {
                    int n_idx = state.idx(ni, nj);
                    double p_n = state.pressures[n_idx];
                    
                    // Upwinding for Oil
                    int up_o = (p_c > p_n) ? c : n_idx;
                    double krog_u, krg_u_o;
                    get_rel_perm(state.gas_saturations[up_o], krog_u, krg_u_o);
                    double To = unit_conv * k_abs * (krog_u / (mu_o * bo)) * (h * (dist == dx ? dy : dx)) / dist;
                    res[2 * c] -= To * (p_n - p_c);

                    // Upwinding for Gas
                    int up_g = (p_c > p_n) ? c : n_idx;
                    double krog_u_g, krg_u;
                    get_rel_perm(state.gas_saturations[up_g], krog_u_g, krg_u);
                    double bg_up = get_bg(state.pressures[up_g]);
                    double Tg = unit_conv * k_abs * (krg_u / (mu_g * bg_up)) * (h * (dist == dx ? dy : dx)) / dist;
                    res[2 * c + 1] -= Tg * (p_n - p_c);
                };

                if (i > 0) add_fluxes(i - 1, j, dx);
                if (i < nx - 1) add_fluxes(i + 1, j, dx);
                if (j > 0) add_fluxes(i, j - 1, dy);
                if (j < ny - 1) add_fluxes(i, j + 1, dy);
            }
        }

        // 3. Source Terms
        for (auto& s : sources) {
            s->apply(res, nullptr, state, dt);
        }

        return res;
    }


    SparseMatrix build_sparse_jacobian(const IState& s_raw, double dt) const override {
        const auto& state = dynamic_cast<const ReservoirOilGas2DState&>(s_raw);
        int nx = state.spatial.nx, ny = state.spatial.ny;
        double dx = state.spatial.dx, dy = state.spatial.dy;
        int n = (int)state.pressures.size();
        
        std::vector<SparseMatrix::Entry> entries;
        double unit_conv = 0.001127;
        double pore_vol_unit = (dx * dy * h) / 5.615;

        #pragma omp parallel
        {
            std::vector<SparseMatrix::Entry> local_entries;
            #pragma omp for
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    int c = state.idx(i, j);
                    double p_c = state.pressures[c];
                    double sg_c = state.gas_saturations[c];
                    double bg_c = get_bg(p_c);
                    double dbg_dp = get_dbg_dp(p_c);

                    // Accumulation
                    local_entries.push_back({2 * c, 2 * c + 1, -phi * pore_vol_unit / (dt * bo)});
                    local_entries.push_back({2 * c + 1, 2 * c, (phi * pore_vol_unit / dt) * sg_c * (-1.0 / (bg_c * bg_c)) * dbg_dp});
                    local_entries.push_back({2 * c + 1, 2 * c + 1, (phi * pore_vol_unit / dt) / bg_c});

                    auto add_flux_jac = [&](int ni, int nj, double dist) {
                        int n_idx = state.idx(ni, nj);
                        double p_n = state.pressures[n_idx];
                        double p_diff = p_n - p_c;
                        
                        int up = (p_c > p_n) ? c : n_idx;
                        double krog, krg;
                        get_rel_perm(state.gas_saturations[up], krog, krg);
                        double bg_up = get_bg(state.pressures[up]);

                        double To = unit_conv * k_abs * (krog / (mu_o * bo)) * (h * (dist == dx ? dy : dx)) / dist;
                        double Tg = unit_conv * k_abs * (krg / (mu_g * bg_up)) * (h * (dist == dx ? dy : dx)) / dist;

                        // Pressure sensitivities
                        local_entries.push_back({2 * c, 2 * c, To});
                        local_entries.push_back({2 * c, 2 * n_idx, -To});
                        local_entries.push_back({2 * c + 1, 2 * c, Tg});
                        local_entries.push_back({2 * c + 1, 2 * n_idx, -Tg});

                        // Saturation sensitivities (Upwind dependency)
                        double sorg_p = 0.1;
                        double s_norm = (1.0 - state.swc - sorg_p);
                        double sge_u = state.gas_saturations[up] / s_norm;
                        
                        double dkrg_dsg = (s_norm > 0 && sge_u > 0 && sge_u < 1) ? (2.0 * sge_u / s_norm) : 0;
                        double dkrog_dsg = (s_norm > 0 && sge_u > 0 && sge_u < 1) ? (-2.0 * (1.0 - sge_u) / s_norm) : 0;

                        double dTo_dsg = unit_conv * k_abs * (dkrog_dsg / (mu_o * bo)) * (h * (dist == dx ? dy : dx)) / dist;
                        double dTg_dsg = unit_conv * k_abs * (dkrg_dsg / (mu_g * bg_up)) * (h * (dist == dx ? dy : dx)) / dist;

                        int up_col = 2 * up + 1;
                        local_entries.push_back({2 * c, up_col, -dTo_dsg * p_diff});
                        local_entries.push_back({2 * c + 1, up_col, -dTg_dsg * p_diff});
                    };

                    if (i > 0) add_flux_jac(i - 1, j, dx);
                    if (i < nx - 1) add_flux_jac(i + 1, j, dx);
                    if (j > 0) add_flux_jac(i, j - 1, dy);
                    if (j < ny - 1) add_flux_jac(i, j + 1, dy);
                }
            }
            #pragma omp critical
            {
                entries.insert(entries.end(), local_entries.begin(), local_entries.end());
            }
        }
        
        Vector dummy_res(2 * n, 0.0);
        for (auto& s : sources) s->apply(dummy_res, nullptr, state, dt, &entries);

        return SparseMatrix::from_triplets(2 * n, 2 * n, entries);
    }


    Matrix build_jacobian(const IState& s_raw, double dt) const override {
        // Fallback for deprecated dense solvers
        SparseMatrix S = build_sparse_jacobian(s_raw, dt);
        Matrix J(2 * S.rows, Vector(2 * S.cols, 0.0)); // Wait, 2*n is already S.rows
        Matrix J_full(S.rows, Vector(S.cols, 0.0));
        for (int i = 0; i < S.rows; ++i) {
            for (int k = S.row_ptr[i]; k < S.row_ptr[i + 1]; ++k) {
                J_full[i][S.col_indices[k]] = S.values[k];
            }
        }
        return J_full;
    }

    Vector evaluate_rhs(const IState& state) const override { return {}; }
};

} // namespace reservoir
} // namespace mod
