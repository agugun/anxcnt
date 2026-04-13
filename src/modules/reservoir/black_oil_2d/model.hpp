#pragma once
#include "lib/modules.hpp"
#include "state.hpp"
#include "../pvt.hpp"
#include "lib/operators.hpp"
#include <vector>

namespace mod {
using namespace top;
namespace reservoir {

/**
 * @brief 2D Black Oil Reservoir Model (3-Phase: Water, Oil, Gas).
 * 
 * Governing Equations:
 *   Conservation of mass for Water, Oil, and Gas phases in 2D space.
 *   Inherits logic from the 3-phase physics but optimized for planar XY flow.
 * 
 * Physical Logic:
 *   - Darcy's Law for multiphase flow.
 *   - Solution gas behavior (Rs) and compressibility handled via BlackOilPVT.
 * 
 * Discretization:
 *   - Spatial: Finite Volume on 5-point Cartesian stencil.
 *   - Upwinding: Single-point upstream weighting for phase mobilities.
 * 
 * Numerical Scheme:
 *   - Fully Implicit (FIM) using Newton-Raphson iteration.
 *   - This module provides robust stability for high-mobility displacements 
 *     or high-contrast PVT behavior compared to IMPES.
 */
class ReservoirBlackOil2DModel : public IModel {
public:
    double k_abs; // md
    double phi;
    double h;     // ft
    
    BlackOilPVT pvt;
    std::vector<std::shared_ptr<mod::ISourceSink>> sources;

    ReservoirBlackOil2DModel(double k, double p, double thickness, 
                             const std::vector<std::shared_ptr<mod::ISourceSink>>& sources_val)
        : k_abs(k), phi(p), h(thickness), sources(sources_val) {}

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
        const auto& state = dynamic_cast<const ReservoirBlackOil2DState&>(s_new);
        const auto& state_old = dynamic_cast<const ReservoirBlackOil2DState&>(s_old);
        int nx = state.spatial.nx, ny = state.spatial.ny;
        double dx = state.spatial.dx, dy = state.spatial.dy;
        int n = nx * ny;
        
        Vector res(3 * n, 0.0);
        double unit_conv = 0.001127;
        double pv_unit = (dx * dy * h) / 5.615;

        #pragma omp parallel for collapse(2)
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                int c = state.idx(i, j);
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

                // 2. Flux
                auto add_flux = [&](int ni, int nj, double dist) {
                    int n_idx = state.idx(ni, nj);
                    double p_n = state.p(n_idx);
                    double area = h * (dist == dx ? dy : dx);
                    double transmissibility = unit_conv * k_abs * area / dist;

                    int up = (p_c > p_n) ? c : n_idx;
                    double krw_u, krog_u, krg_u;
                    get_rel_perm(state.sw(up), state.sg(up), krw_u, krog_u, krg_u);

                    double p_up = state.p(up);
                    double rs_up = pvt.get_rs(p_up);
                    double bo_up = pvt.get_bo(p_up, rs_up);
                    double bg_up = pvt.get_bg(p_up);
                    double bw_up = pvt.get_bw(p_up);

                    res[3 * c]     -= transmissibility * (krw_u / (pvt.get_mu_w(p_up) * bw_up)) * (p_n - p_c);
                    res[3 * c + 1] -= transmissibility * (krog_u / (pvt.get_mu_o(p_up, rs_up) * bo_up)) * (p_n - p_c);
                    double lam_g = krg_u / (pvt.get_mu_g(p_up) * bg_up);
                    double lam_rs = rs_up * krog_u / (pvt.get_mu_o(p_up, rs_up) * bo_up);
                    res[3 * c + 2] -= transmissibility * (lam_g + lam_rs) * (p_n - p_c);
                };

                if (i > 0) add_flux(i - 1, j, dx);
                if (i < nx - 1) add_flux(i + 1, j, dx);
                if (j > 0) add_flux(i, j - 1, dy);
                if (j < ny - 1) add_flux(i, j + 1, dy);
            }
        }
        for (auto& s : sources) s->apply(res, nullptr, state, dt);
        return res;
    }

    SparseMatrix build_sparse_jacobian(const IState& s_raw, double dt) const override {
        const auto& state = dynamic_cast<const ReservoirBlackOil2DState&>(s_raw);
        int nx = state.spatial.nx, ny = state.spatial.ny;
        double dx = state.spatial.dx, dy = state.spatial.dy;
        int n = nx * ny;
        
        std::vector<SparseMatrix::Entry> entries;
        double pv_unit = (phi * dx * dy * h) / 5.615;
        double unit_conv = 0.001127;

        #pragma omp parallel
        {
            std::vector<SparseMatrix::Entry> local_entries;
            #pragma omp for
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    int c = state.idx(i, j);
                    double p = state.p(c), sw = state.sw(c), sg = state.sg(c), so = state.so(c);
                    double rs = pvt.get_rs(p), bo = pvt.get_bo(p, rs), bg = pvt.get_bg(p), bw = pvt.get_bw(p);
                    
                    // --- Accumulation ---
                    double cw = 3e-6; 
                    double d1bw_dp = cw / bw;
                    double d1bo_dp = -1.0/(bo*bo) * pvt.get_dbo_dp(p);
                    double d1bg_dp = -1.0/(bg*bg) * pvt.get_dbg_dp(p);
                    double drs_dp  = pvt.get_drs_dp(p);

                    local_entries.push_back({3 * c, 3 * c, (phi*pv_unit/dt) * sw * d1bw_dp});
                    local_entries.push_back({3 * c, 3 * c + 1, (phi*pv_unit/dt) / bw});
                    
                    local_entries.push_back({3 * c + 1, 3 * c, (phi*pv_unit/dt) * so * d1bo_dp});
                    local_entries.push_back({3 * c + 1, 3 * c + 1, (phi*pv_unit/dt) * (-1.0/bo)});
                    local_entries.push_back({3 * c + 1, 3 * c + 2, (phi*pv_unit/dt) * (-1.0/bo)});

                    double dAccG_dp = sg*d1bg_dp + so*(drs_dp/bo + rs*d1bo_dp);
                    local_entries.push_back({3 * c + 2, 3 * c, (phi*pv_unit/dt) * dAccG_dp});
                    local_entries.push_back({3 * c + 2, 3 * c + 1, (phi*pv_unit/dt) * (-rs/bo)});
                    local_entries.push_back({3 * c + 2, 3 * c + 2, (phi*pv_unit/dt) * (1.0/bg - rs/bo)});

                    // --- Flux ---
                    auto add_flux_jac = [&](int ni, int nj, double dist) {
                        int n_idx = state.idx(ni, nj);
                        double p_n = state.p(n_idx);
                        double p_diff = p_n - p;
                        double area = h * (dist == dx ? dy : dx);
                        double trans = unit_conv * k_abs * area / dist;

                        int up = (p > p_n) ? c : n_idx;
                        double krw, krog, krg;
                        get_rel_perm(state.sw(up), state.sg(up), krw, krog, krg);
                        double p_up = state.p(up);
                        double rs_up = pvt.get_rs(p_up);

                        double muw_u = pvt.get_mu_w(p_up), bw_u = pvt.get_bw(p_up);
                        double muo_u = pvt.get_mu_o(p_up, rs_up), bo_u = pvt.get_bo(p_up, rs_up);
                        double mug_u = pvt.get_mu_g(p_up), bg_u = pvt.get_bg(p_up);

                        double Tw = trans * krw / (muw_u * bw_u);
                        double To = trans * krog / (muo_u * bo_u);
                        double Tg = trans * (krg / (mug_u * bg_u) + rs_up * krog / (muo_u * bo_u));

                        local_entries.push_back({3 * c, 3 * c, Tw});
                        local_entries.push_back({3 * c, 3 * n_idx, -Tw});
                        local_entries.push_back({3 * c + 1, 3 * c, To});
                        local_entries.push_back({3 * c + 1, 3 * n_idx, -To});
                        local_entries.push_back({3 * c + 2, 3 * c, Tg});
                        local_entries.push_back({3 * c + 2, 3 * n_idx, -Tg});

                        // Saturation sensitivities
                        double swc_p = 0.2, sorg_p = 0.1, sgr_p = 0.05, swr_p = 0.2;
                        double sw_u = state.sw(up), sg_u = state.sg(up), so_u = 1.0 - sw_u - sg_u;
                        double swe = (sw_u - swr_p) / (1.0 - swr_p - sorg_p);
                        double sge = (sg_u - sgr_p) / (1.0 - swc_p - sgr_p);
                        double soe = (so_u - sorg_p) / (1.0 - swr_p - sorg_p);

                        double dkrw_dsw = (swe > 0 && swe < 1) ? 2.0 * swe / (1.0 - swr_p - sorg_p) : 0;
                        double dkrg_dsg = (sge > 0 && sge < 1) ? 2.0 * sge / (1.0 - swc_p - sgr_p) : 0;
                        double dkro_dso = (soe > 0 && soe < 1) ? 3.0 * soe * soe / (1.0 - swr_p - sorg_p) : 0;

                        double dTw_dsw = trans * (dkrw_dsw / (muw_u * bw_u));
                        double dTo_dsw = trans * (-dkro_dso / (muo_u * bo_u));
                        double dTo_dsg = trans * (-dkro_dso / (muo_u * bo_u));
                        double dTg_dsw = trans * (rs_up * -dkro_dso / (muo_u * bo_u));
                        double dTg_dsg = trans * (dkrg_dsg / (mug_u * bg_u) + rs_up * -dkro_dso / (muo_u * bo_u));

                        int up_col_sw = 3 * up + 1, up_col_sg = 3 * up + 2;
                        local_entries.push_back({3 * c, up_col_sw, -dTw_dsw * p_diff});
                        local_entries.push_back({3 * c + 1, up_col_sw, -dTo_dsw * p_diff});
                        local_entries.push_back({3 * c + 1, up_col_sg, -dTo_dsg * p_diff});
                        local_entries.push_back({3 * c + 2, up_col_sw, -dTg_dsw * p_diff});
                        local_entries.push_back({3 * c + 2, up_col_sg, -dTg_dsg * p_diff});
                    };
                    if (i > 0)      add_flux_jac(i - 1, j, dx);
                    if (i < nx - 1) add_flux_jac(i + 1, j, dx);
                    if (j > 0)      add_flux_jac(i, j - 1, dy);
                    if (j < ny - 1) add_flux_jac(i, j + 1, dy);
                }
            }
            #pragma omp critical
            {
                entries.insert(entries.end(), local_entries.begin(), local_entries.end());
            }
        }
        Vector dummy_res(3 * n, 0.0);
        for (auto& s : sources) s->apply(dummy_res, nullptr, state, dt, &entries);
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
    Vector evaluate_rhs(const IState& state) const override { return {}; }
};
} // namespace reservoir
} // namespace mod
