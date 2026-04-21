#pragma once
#include "lib/interfaces.hpp"
#include "state.hpp"
#include "lib/discretization.hpp"
#include <vector>

namespace mod::reservoir {
using namespace top;

/**
 * @brief 2D Oil-Gas Reservoir Physical Model (Properties).
 */
class OilGas2DModel : public IModel {
public:
    double mu_o, mu_g; 
    double bo;         // Oil formation volume factor (assumed constant)
    double p_std, bg_std; 
    std::shared_ptr<num::discretization::Conductance2D> rock_cond;
    double pore_vol_per_cell; 
    std::vector<std::shared_ptr<ISourceSink>> wells;

    OilGas2DModel(std::shared_ptr<num::discretization::Conductance2D> cond, double pv, 
                  double mo, double mg, double b_o, double ps, double bgs,
                  const std::vector<std::shared_ptr<ISourceSink>>& wells_val)
        : rock_cond(cond), pore_vol_per_cell(pv), mu_o(mo), mu_g(mg), 
          bo(b_o), p_std(ps), bg_std(bgs), wells(wells_val) {}

    double get_bg(double p) const { return bg_std * (p_std / std::max(1.0, p)); }
    double get_tolerance() const override { return 1e-4; }

    void get_rel_perm(double sg, double& krog, double& krg) const {
        double swc = 0.2, sorg = 0.1;
        double sge = sg / (1.0 - swc - sorg);
        sge = std::max(0.0, std::min(1.0, sge));
        krg = sge * sge;
        krog = (1.0 - sge) * (1.0 - sge);
    }

    Vector get_accumulation_weights(const IGrid& grid, const IState& state) const override {
        size_t n = state.to_vector().size();
        return Vector(n, pore_vol_per_cell);
    }
};

/**
 * @brief 2D Oil-Gas FVM Discretizer (Numerical Assembly).
 */
class OilGas2DDiscretizer : public IDiscretizer {
public:
    void assemble_jacobian(const IGrid& grid, const IModel& model, const IState& state, SparseMatrix& J) const override {
        const auto& m = static_cast<const OilGas2DModel&>(model);
        const auto& s = static_cast<const ReservoirOilGas2DState&>(state);
        int nx = (int)s.spatial->nx;
        int ny = (int)s.spatial->ny;
        int n = nx * ny;

        if (J.rows != 2 * n) J = SparseMatrix(2 * n, 2 * n);
        J.triplets.clear();

        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                int c = s.spatial->idx(i, j);
                double p = s.p(c), sg = s.sg(c);
                double bg = m.get_bg(p);

                // Diagonal accumulation block sensitivities
                J.triplets.push_back({2 * c, 2 * c + 1, -1.0 / m.bo}); // dRo/dSg (approx)
                J.triplets.push_back({2 * c + 1, 2 * c + 1, 1.0 / bg}); // dRg/dSg

                auto add_flux_jac = [&](int n_idx, double t_rock) {
                    int up = (s.p(c) > s.p(n_idx)) ? c : n_idx;
                    double krog, krg;
                    m.get_rel_perm(s.sg(up), krog, krg);
                    double bgu = m.get_bg(s.p(up));

                    double To = t_rock * (krog / (m.mu_o * m.bo));
                    double Tg = t_rock * (krg / (m.mu_g * bgu));

                    J.triplets.push_back({2 * c, 2 * c, To});
                    J.triplets.push_back({2 * c, 2 * n_idx, -To});
                    J.triplets.push_back({2 * c + 1, 2 * c, Tg});
                    J.triplets.push_back({2 * c + 1, 2 * n_idx, -Tg});
                };

                if (i > 0)      add_flux_jac(s.spatial->idx(i - 1, j), m.rock_cond->Tx[j * (nx - 1) + i - 1]);
                if (i < nx - 1) add_flux_jac(s.spatial->idx(i + 1, j), m.rock_cond->Tx[j * (nx - 1) + i]);
                if (j > 0)      add_flux_jac(s.spatial->idx(i, j - 1), m.rock_cond->Ty[(j - 1) * nx + i]);
                if (j < ny - 1) add_flux_jac(s.spatial->idx(i, j + 1), m.rock_cond->Ty[j * nx + i]);
            }
        }
    }

    void assemble_residual(const IGrid& grid, const IModel& model, const IState& state, Vector& R) const override {
        const auto& m = static_cast<const OilGas2DModel&>(model);
        const auto& s = static_cast<const ReservoirOilGas2DState&>(state);
        int nx = (int)s.spatial->nx;
        int ny = (int)s.spatial->ny;

        #pragma omp parallel for collapse(2)
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                int c = s.spatial->idx(i, j);
                double p_c = s.p(c);
                double net_o = 0, net_g = 0;

                auto add_flux_res = [&](int n_idx, double t_rock) {
                    double p_n = s.p(n_idx);
                    int up = (p_c > p_n) ? c : n_idx;
                    double krog_u, krg_u;
                    m.get_rel_perm(s.sg(up), krog_u, krg_u);
                    double bgu = m.get_bg(s.p(up));

                    net_o += t_rock * (krog_u / (m.mu_o * m.bo)) * (p_n - p_c);
                    net_g += t_rock * (krg_u / (m.mu_g * bgu)) * (p_n - p_c);
                };

                if (i > 0)      add_flux_res(s.spatial->idx(i - 1, j), m.rock_cond->Tx[j * (nx - 1) + i - 1]);
                if (i < nx - 1) add_flux_res(s.spatial->idx(i + 1, j), m.rock_cond->Tx[j * (nx - 1) + i]);
                if (j > 0)      add_flux_res(s.spatial->idx(i, j - 1), m.rock_cond->Ty[(j - 1) * nx + i]);
                if (j < ny - 1) add_flux_res(s.spatial->idx(i, j + 1), m.rock_cond->Ty[j * nx + i]);

                R[2*c] = -net_o;
                R[2*c+1] = -net_g;
            }
        }
    }

    void apply_boundary_conditions(const IGrid& grid, const IModel& model, const IState& state, SparseMatrix& J, Vector& R) const override {
        // No-flow is default.
    }
};

} // namespace mod::reservoir
