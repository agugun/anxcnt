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

        #pragma omp parallel for collapse(3)
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
                        double p_diff = p_n - p_c;
                        
                        int up = (p_c > p_n) ? c : n_idx; 
                        double krw_u, krog_u, krg_u;
                        get_rel_perm(state.sw(up), state.sg(up), krw_u, krog_u, krg_u);

                        double p_up = state.p(up);
                        double rs_up = pvt.get_rs(p_up);
                        double trans = unit_conv * k_abs * area / dist;

                        res[3 * c]     -= trans * (krw_u / (pvt.get_mu_w(p_up) * pvt.get_bw(p_up))) * p_diff;
                        res[3 * c + 1] -= trans * (krog_u / (pvt.get_mu_o(p_up, rs_up) * pvt.get_bo(p_up, rs_up))) * p_diff;
                        double lam_g = krg_u / (pvt.get_mu_g(p_up) * pvt.get_bg(p_up));
                        double lam_rs = rs_up * krog_u / (pvt.get_mu_o(p_up, rs_up) * pvt.get_bo(p_up, rs_up));
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
        #pragma omp parallel
        {
            std::vector<SparseMatrix::Entry> local_entries;
            #pragma omp for
            for (int k = 0; k < nz; ++k) {
                for (int j = 0; j < ny; ++j) {
                    for (int i = 0; i < nx; ++i) {
                        int c = state.idx(i, j, k);
                        double p_c = state.p(c);
                        double sw_c = state.sw(c), sg_c = state.sg(c);
                        
                        // Local Accumulation Sensitivities
                        double bw_c = pvt.get_bw(p_c), bo_c = pvt.get_bo(p_c, pvt.get_rs(p_c)), bg_c = pvt.get_bg(p_c);
                        double dbw_dp = pvt.get_dbw_dp(p_c), dbo_dp = pvt.get_dbo_dp(p_c), dbg_dp = pvt.get_dbg_dp(p_c);
                        double rs_c = pvt.get_rs(p_c), drs_dp = pvt.get_drs_dp(p_c);
                        double pv_term = phi * (state.spatial.dx * state.spatial.dy * state.spatial.dz) / (5.615 * dt);

                        // Mass balance derivatives
                        local_entries.push_back(SparseMatrix::Entry{3 * c, 3 * c, pv_term * sw_c * (-1.0 / (bw_c * bw_c)) * dbw_dp});
                        local_entries.push_back(SparseMatrix::Entry{3 * c, 3 * c + 1, pv_term / bw_c});
                        
                        local_entries.push_back(SparseMatrix::Entry{3 * c + 1, 3 * c, pv_term * (1.0 - sw_c - sg_c) * (-1.0 / (bo_c * bo_c)) * dbo_dp});
                        local_entries.push_back(SparseMatrix::Entry{3 * c + 1, 3 * c + 1, -pv_term / bo_c});
                        local_entries.push_back(SparseMatrix::Entry{3 * c + 1, 3 * c + 2, -pv_term / bo_c});

                        double mass_g_dp = pv_term * (sg_c * (-1.0 / (bg_c * bg_c)) * dbg_dp + rs_c * (1.0 - sw_c - sg_c) * (-1.0 / (bo_c * bo_c)) * dbo_dp + drs_dp * (1.0 - sw_c - sg_c) / bo_c);
                        local_entries.push_back(SparseMatrix::Entry{3 * c + 2, 3 * c, mass_g_dp});
                        local_entries.push_back(SparseMatrix::Entry{3 * c + 2, 3 * c + 1, -pv_term * rs_c / bo_c});
                        local_entries.push_back(SparseMatrix::Entry{3 * c + 2, 3 * c + 2, pv_term * (1.0 / bg_c - rs_c / bo_c)});

                        auto add_flux_jac = [&](int ni, int nj, int nk, double dist, double area) {
                            int n_idx = state.idx(ni, nj, nk);
                            double p_n = state.p(n_idx);
                            double p_diff = p_n - p_c;
                            int up = (p_c > p_n) ? c : n_idx;

                            double p_up = state.p(up);
                            double muw_u = pvt.get_mu_w(p_up), muo_u = pvt.get_mu_o(p_up, pvt.get_rs(p_up)), mug_u = pvt.get_mu_g(p_up);
                            double bw_u = pvt.get_bw(p_up), bo_u = pvt.get_bo(p_up, pvt.get_rs(p_up)), bg_u = pvt.get_bg(p_up);
                            double rs_u = pvt.get_rs(p_up);
                            
                            double krw, krog, krg;
                            get_rel_perm(state.sw(up), state.sg(up), krw, krog, krg);

                            double trans = 0.001127 * k_abs * area / dist;
                            double Tw = trans * (krw / (muw_u * bw_u));
                            double To = trans * (krog / (muo_u * bo_u));
                            double Tg = trans * (krg / (mug_u * bg_u) + rs_u * krog / (muo_u * bo_u));

                            local_entries.push_back(SparseMatrix::Entry{3 * c, 3 * c, Tw});
                            local_entries.push_back(SparseMatrix::Entry{3 * c, 3 * n_idx, -Tw});
                            local_entries.push_back(SparseMatrix::Entry{3 * c + 1, 3 * c, To});
                            local_entries.push_back(SparseMatrix::Entry{3 * c + 1, 3 * n_idx, -To});
                            local_entries.push_back(SparseMatrix::Entry{3 * c + 2, 3 * c, Tg});
                            local_entries.push_back(SparseMatrix::Entry{3 * c + 2, 3 * n_idx, -Tg});

                            // Saturation derivatives (Upwind dependency)
                            double swc_p = 0.2, sorg_p = 0.1, sgr_p = 0.05, swr_p = 0.2;
                            double sw_u = state.sw(up), sg_u = state.sg(up);
                            double so_u = 1.0 - sw_u - sg_u;
                            
                            double w_norm = (1.0 - swr_p - sorg_p);
                            double g_norm = (1.0 - swc_p - sgr_p);
                            double o_norm = (1.0 - swr_p - sorg_p);

                            double swe_curr = (sw_u - swr_p) / w_norm;
                            double sge_u = (sg_u - sgr_p) / g_norm;
                            double soe_u = (so_u - sorg_p) / o_norm;

                            double dkrw_dsw = (w_norm > 0 && swe_curr > 0 && swe_curr < 1) ? (2.0 * swe_curr / w_norm) : 0;
                            double dkrg_dsg = (g_norm > 0 && sge_u > 0 && sge_u < 1) ? (2.0 * sge_u / g_norm) : 0;
                            double dkro_dso = (o_norm > 0 && soe_u > 0 && soe_u < 1) ? (3.0 * soe_u * soe_u / o_norm) : 0;

                            double dTw_dsw = trans * (dkrw_dsw / (muw_u * bw_u));
                            double dTo_dsw = trans * (dkro_dso * -1.0 / (muo_u * bo_u));
                            double dTo_dsg = trans * (dkro_dso * -1.0 / (muo_u * bo_u));
                            double dTg_dsw = trans * (rs_u * dkro_dso * -1.0 / (muo_u * bo_u));
                            double dTg_dsg = trans * (dkrg_dsg / (mug_u * bg_u) + rs_u * dkro_dso * -1.0 / (muo_u * bo_u));

                            int col_sw = 3 * up + 1, col_sg = 3 * up + 2;
                            local_entries.push_back(SparseMatrix::Entry{3 * c, col_sw, -dTw_dsw * p_diff});
                            local_entries.push_back(SparseMatrix::Entry{3 * c + 1, col_sw, -dTo_dsw * p_diff});
                            local_entries.push_back(SparseMatrix::Entry{3 * c + 1, col_sg, -dTo_dsg * p_diff});
                            local_entries.push_back(SparseMatrix::Entry{3 * c + 2, col_sw, -dTg_dsw * p_diff});
                            local_entries.push_back(SparseMatrix::Entry{3 * c + 2, col_sg, -dTg_dsg * p_diff});
                        };

                        if (i > 0)      add_flux_jac(i - 1, j, k, dx, dy * dz);
                        if (i < nx - 1) add_flux_jac(i + 1, j, k, dx, dy * dz);
                        if (j > 0)      add_flux_jac(i, j - 1, k, dy, dx * dz);
                        if (j < ny - 1) add_flux_jac(i, j + 1, k, dy, dx * dz);
                        if (k > 0)      add_flux_jac(i, j, k - 1, dz, dx * dy);
                        if (k < nz - 1) add_flux_jac(i, j, k + 1, dz, dx * dy);
                    }
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
