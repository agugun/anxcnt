#pragma once
#include "lib/modules.hpp"
#include "state.hpp"
#include "../pvt.hpp"
#include "lib/operators.hpp"
#include <vector>
#include <cmath>

namespace mod {
using namespace top;
namespace reservoir {

/**
 * @brief 3D Black Oil Reservoir Model (3-Phase: Water, Oil, Gas).
 * 
 * Governing Equations (Mass Conservation):
 *   1. Water: div( (k*krw)/(muw*Bw) * grad(P) ) = d(phi*Sw/Bw)/dt
 *   2. Oil:   div( (k*kro)/(muo*Bo) * grad(P) ) = d(phi*So/Bo)/dt
 *   3. Gas:   div( (k*krg)/(mug*Bg) * grad(P) + Rs*(k*kro)/(muo*Bo) * grad(P) ) 
 *             = d(phi*(Sg/Bg + Rs*So/Bo))/dt
 * 
 * Physical Logic:
 *   - Multiphase Flow: Darcy's Law extended with relative permeabilities.
 *   - Solution Gas: Oil can contain dissolved gas (Rs), which releases as 
 *     pressure drops below bubble point.
 *   - Gravity: (Currently simplified or neglected in this stencil).
 * 
 * Discretization:
 *   - Spatial: Finite Volume Method (FVM) with a 7-point orthogonal stencil (3D).
 *   - Upwinding: Phase mobilities are evaluated using Single-Point Upwind 
 *     based on potential (pressure) gradients between cells.
 * 
 * Numerical Scheme:
 *   - Fully Implicit (FIM): Both Pressure and Saturations are solved simultaneously.
 *   - Newton-Raphson: Used to resolve the non-linear residuals.
 * 
 * Jacobian:
 *   Constructed with 3x3 blocks per cell. Includes sensitivities for:
 *   - Accumulation: Volume changes due to compressibility and saturation.
 *   - Flux: Mobility and potential changes between neighboring cells.
 *   - Coupling: Strong cross-phase coupling between P, Sw, and Sg.
 */
class ReservoirBlackOil3DModel : public IModel {
public:
    double k_abs; // md
    double phi;
    
    BlackOilPVT pvt;
    std::vector<std::shared_ptr<mod::ISourceSink>> sources;

    ReservoirBlackOil3DModel(double k, double p, 
                             const std::vector<std::shared_ptr<mod::ISourceSink>>& sources_val)
        : k_abs(k), phi(p), sources(sources_val) {}

    // 3-Phase Relative Permeability (Stone's Method II Proxy)
    void get_rel_perm(double sw, double sg, double& krw, double& krog, double& krg) const {
        double swc = 0.2, sorg = 0.1, sgr = 0.05, swr = 0.2;
        
        // Water
        double swe = (sw - swr) / (1.0 - swr - sorg);
        krw = std::pow(std::max(0.0, std::min(1.0, swe)), 2.0);
        
        // Gas
        double sge = (sg - sgr) / (1.0 - swc - sgr);
        krg = std::pow(std::max(0.0, std::min(1.0, sge)), 2.0);
        
        // Oil (Stone II simplified)
        double so = 1.0 - sw - sg;
        double soe = (so - sorg) / (1.0 - swr - sorg);
        krog = std::pow(std::max(0.0, std::min(1.0, soe)), 3.0);
    }

    Vector build_residual(const IState& s_new, const IState& s_old, double dt) const override {
        const auto& state = dynamic_cast<const ReservoirBlackOil3DState&>(s_new);
        const auto& state_old = dynamic_cast<const ReservoirBlackOil3DState&>(s_old);
        int nx = state.spatial.nx, ny = state.spatial.ny, nz = state.spatial.nz;
        double dx = state.spatial.dx, dy = state.spatial.dy, dz = state.spatial.dz;
        int n = nx * ny * nz;
        
        Vector res(3 * n, 0.0);
        double unit_conv = 0.001127;
        double cell_vol = dx * dy * dz;
        double pv_unit = cell_vol / 5.615;

        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    int c = state.idx(i, j, k);
                    double p_c = state.p(c);
                    double sw_c = state.sw(c), sg_c = state.sg(c), so_c = state.so(c);
                    
                    double rs_c = pvt.get_rs(p_c);
                    double bo_c = pvt.get_bo(p_c, rs_c);
                    double bg_c = pvt.get_bg(p_c);
                    double bw_c = pvt.get_bw(p_c);

                    // 1. Accumulation
                    res[3 * c]     = (phi * pv_unit / dt) * (sw_c / bw_c - state_old.sw(c) / pvt.get_bw(state_old.p(c)));
                    res[3 * c + 1] = (phi * pv_unit / dt) * (so_c / bo_c - state_old.so(c) / pvt.get_bo(state_old.p(c), pvt.get_rs(state_old.p(c))));
                    double acc_g = sg_c / bg_c + rs_c * so_c / bo_c;
                    double acc_g_old = state_old.sg(c) / pvt.get_bg(state_old.p(c)) + pvt.get_rs(state_old.p(c)) * state_old.so(c) / pvt.get_bo(state_old.p(c), pvt.get_rs(state_old.p(c)));
                    res[3 * c + 2] = (phi * pv_unit / dt) * (acc_g - acc_g_old);

                    // 2. Flux (6 directions)
                    auto add_flux = [&](int ni, int nj, int nk, double area, double dist, double dz_grad) {
                        int n_idx = state.idx(ni, nj, nk);
                        double p_n = state.p(n_idx);
                        
                        // Gravity head (simplified: water=0.433 psi/ft, oil=0.35, gas=0.08)
                        // Pot_diff = (P_n - P_c) - rho_g * (z_n - z_c)
                        // Here dz_grad is (nk_n - nk_c) * dz
                        double p_diff = p_n - p_c;
                        
                        int up = (p_c > p_n) ? c : n_idx; // Pressure based upwinding for simplicity
                        double krw, krog, krg;
                        get_rel_perm(state.sw(up), state.sg(up), krw, krog, krg);

                        double p_up = state.p(up);
                        double rs_up = pvt.get_rs(p_up);
                        double trans = unit_conv * k_abs * area / dist;

                        res[3 * c]     -= trans * (krw / (pvt.get_mu_w(p_up) * pvt.get_bw(p_up))) * p_diff;
                        res[3 * c + 1] -= trans * (krog / (pvt.get_mu_o(p_up, rs_up) * pvt.get_bo(p_up, rs_up))) * p_diff;
                        double lam_g = krg / (pvt.get_mu_g(p_up) * pvt.get_bg(p_up));
                        double lam_rs = rs_up * krog / (pvt.get_mu_o(p_up, rs_up) * pvt.get_bo(p_up, rs_up));
                        res[3 * c + 2] -= trans * (lam_g + lam_rs) * p_diff;
                    };

                    if (i > 0)      add_flux(i - 1, j, k, dy * dz, dx, 0.0);
                    if (i < nx - 1) add_flux(i + 1, j, k, dy * dz, dx, 0.0);
                    if (j > 0)      add_flux(i, j - 1, k, dx * dz, dy, 0.0);
                    if (j < ny - 1) add_flux(i, j + 1, k, dx * dz, dy, 0.0);
                    if (k > 0)      add_flux(i, j, k - 1, dx * dy, dz, -dz);
                    if (k < nz - 1) add_flux(i, j, k + 1, dx * dy, dz, dz);
                }
            }
        }
        for (auto& s : sources) s->apply(res, nullptr, state, dt);
        return res;
    }

    SparseMatrix build_sparse_jacobian(const IState& s_raw, double dt) const override {
        const auto& state = dynamic_cast<const ReservoirBlackOil3DState&>(s_raw);
        int nx = state.spatial.nx, ny = state.spatial.ny, nz = state.spatial.nz;
        double dx = state.spatial.dx, dy = state.spatial.dy, dz = state.spatial.dz;
        int n = nx * ny * nz;

        std::vector<SparseMatrix::Entry> entries;
        double pv_unit = (phi * dx * dy * dz) / 5.615;
        double unit_conv = 0.001127;

        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    int c = state.idx(i, j, k);
                    double p = state.p(c), sw = state.sw(c), sg = state.sg(c), so = state.so(c);
                    double rs = pvt.get_rs(p), bo = pvt.get_bo(p, rs), bg = pvt.get_bg(p), bw = pvt.get_bw(p);
                    
                    // --- Accumulation Diagonals ---
                    double cw = 3e-6; 
                    double d1bw_dp = cw / bw;
                    double d1bo_dp = -1.0/(bo*bo) * pvt.get_dbo_dp(p);
                    double d1bg_dp = -1.0/(bg*bg) * pvt.get_dbg_dp(p);
                    double drs_dp  = pvt.get_drs_dp(p);

                    entries.push_back({3 * c, 3 * c, (phi*pv_unit/dt) * sw * d1bw_dp});
                    entries.push_back({3 * c, 3 * c + 1, (phi*pv_unit/dt) / bw});
                    
                    entries.push_back({3 * c + 1, 3 * c, (phi*pv_unit/dt) * so * d1bo_dp});
                    entries.push_back({3 * c + 1, 3 * c + 1, (phi*pv_unit/dt) * (-1.0/bo)});
                    entries.push_back({3 * c + 1, 3 * c + 2, (phi*pv_unit/dt) * (-1.0/bo)});

                    double dAccG_dp = sg*d1bg_dp + so*(drs_dp/bo + rs*d1bo_dp);
                    entries.push_back({3 * c + 2, 3 * c, (phi*pv_unit/dt) * dAccG_dp});
                    entries.push_back({3 * c + 2, 3 * c + 1, (phi*pv_unit/dt) * (-rs/bo)});
                    entries.push_back({3 * c + 2, 3 * c + 2, (phi*pv_unit/dt) * (1.0/bg - rs/bo)});

                    // --- Flux Neighbors ---
                    auto add_flux_jac = [&](int ni, int nj, int nk, double area, double dist) {
                        int n_idx = state.idx(ni, nj, nk);
                        double p_n = state.p(n_idx);
                        double p_diff = p_n - p;
                        double trans = unit_conv * k_abs * area / dist;

                        int up = (p > p_n) ? c : n_idx;
                        double krw, krog, krg;
                        get_rel_perm(state.sw(up), state.sg(up), krw, krog, krg);
                        double p_u = state.p(up);
                        double rs_u = pvt.get_rs(p_u);

                        double muw_u = pvt.get_mu_w(p_u), bw_u = pvt.get_bw(p_u);
                        double muo_u = pvt.get_mu_o(p_u, rs_u), bo_u = pvt.get_bo(p_u, rs_u);
                        double mug_u = pvt.get_mu_g(p_u), bg_u = pvt.get_bg(p_u);

                        double Tw = trans * krw / (muw_u * bw_u);
                        double To = trans * krog / (muo_u * bo_u);
                        double Tg = trans * (krg / (mug_u * bg_u) + rs_u * krog / (muo_u * bo_u));

                        // Pressure sensitivities (Transmissibility * dP)
                        entries.push_back({3 * c, 3 * c, Tw});
                        entries.push_back({3 * c, 3 * n_idx, -Tw});
                        entries.push_back({3 * c + 1, 3 * c, To});
                        entries.push_back({3 * c + 1, 3 * n_idx, -To});
                        entries.push_back({3 * c + 2, 3 * c, Tg});
                        entries.push_back({3 * c + 2, 3 * n_idx, -Tg});

                        // Saturation derivatives (Upwind dependency)
                        double swc_p = 0.2, sorg_p = 0.1, sgr_p = 0.05, swr_p = 0.2;
                        double sw_u = state.sw(up), sg_u = state.sg(up);
                        double so_u = 1.0 - sw_u - sg_u;
                        
                        double w_norm = (1.0 - swr_p - sorg_p);
                        double g_norm = (1.0 - swc_p - sgr_p);
                        double o_norm = (1.0 - swr_p - sorg_p);

                        double swe_curr = (sw_u - swr_p) / w_norm;
                        double sge_curr = (sg_u - sgr_p) / g_norm;
                        double soe_curr = (so_u - sorg_p) / o_norm;

                        double dkrw_dsw = (w_norm > 0 && swe_curr > 0 && swe_curr < 1) ? (2.0 * swe_curr / w_norm) : 0;
                        double dkrg_dsg = (g_norm > 0 && sge_curr > 0 && sge_curr < 1) ? (2.0 * sge_curr / g_norm) : 0;
                        double dkro_dso = (o_norm > 0 && soe_curr > 0 && soe_curr < 1) ? (3.0 * soe_curr * soe_curr / o_norm) : 0;

                        double dTw_dsw = trans * (dkrw_dsw / (muw_u * bw_u));
                        double dTo_dsw = trans * (dkro_dso * -1.0 / (muo_u * bo_u)); // dTo/dSo * dSo/dSw
                        double dTo_dsg = trans * (dkro_dso * -1.0 / (muo_u * bo_u)); // dTo/dSo * dSo/dSg
                        
                        double dTg_dsw = trans * (rs_u * dkro_dso * -1.0 / (muo_u * bo_u));
                        double dTg_dsg = trans * (dkrg_dsg / (mug_u * bg_u) + rs_u * dkro_dso * -1.0 / (muo_u * bo_u));

                        // Add to correct columns based on upwind cell
                        int col_sw = 3 * up + 1;
                        int col_sg = 3 * up + 2;
                        entries.push_back({3 * c, col_sw, -dTw_dsw * p_diff});
                        entries.push_back({3 * c + 1, col_sw, -dTo_dsw * p_diff});
                        entries.push_back({3 * c + 1, col_sg, -dTo_dsg * p_diff});
                        entries.push_back({3 * c + 2, col_sw, -dTg_dsw * p_diff});
                        entries.push_back({3 * c + 2, col_sg, -dTg_dsg * p_diff});
                    };
                    if (i > 0)      add_flux_jac(i - 1, j, k, dy * dz, dx);
                    if (i < nx - 1) add_flux_jac(i + 1, j, k, dy * dz, dx);
                    if (j > 0)      add_flux_jac(i, j - 1, k, dx * dz, dy);
                    if (j < ny - 1) add_flux_jac(i, j + 1, k, dx * dz, dy);
                    if (k > 0)      add_flux_jac(i, j, k - 1, dx * dy, dz);
                    if (k < nz - 1) add_flux_jac(i, j, k + 1, dx * dy, dz);
                }
            }
        }
        Vector dummy_res(3 * n, 0.0);
        Matrix dummy_jac(3 * n, Vector(3 * n, 0.0));
        for (auto& s : sources) s->apply(dummy_res, &dummy_jac, state, dt);
        for(int r=0; r<3*n; ++r) {
            for(int cl=0; cl<3*n; ++cl) {
                if(dummy_jac[r][cl] != 0) entries.push_back({r, cl, dummy_jac[r][cl]});
            }
        }
        return SparseMatrix::from_triplets(3 * n, 3 * n, entries);
    }

    Matrix build_jacobian(const IState& s_raw, double dt) const override {
        SparseMatrix S = build_sparse_jacobian(s_raw, dt);
        Matrix J(S.rows, Vector(S.cols, 0.0));
        for (int i = 0; i < S.rows; ++i) {
            for (int k = S.row_ptr[i]; k < S.row_ptr[i + 1]; ++k) {
                J[i][S.col_indices[k]] = S.values[k];
            }
        }
        return J;
    }
    Vector apply_jacobian(const IState& state_raw, const Vector& v, double dt) const override {
        Matrix J = build_jacobian(state_raw, dt);
        int rows = (int)J.size();
        int cols = (int)J[0].size();
        Vector product(rows, 0.0);
        for (int i = 0; i < rows; ++i) {
            for (int k = 0; k < cols; ++k) {
                product[i] += J[i][k] * v[k];
            }
        }
        return product;
    }

    Vector evaluate_rhs(const IState& state) const override { return {}; }
};

} // namespace reservoir
} // namespace mod
