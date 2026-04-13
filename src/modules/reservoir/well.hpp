#pragma once
#include "lib/operators.hpp"
#include <vector>
#include <memory>
#include <cmath>
#include <functional>

namespace mod {
using namespace top;

/**
 * @brief Base class for Grid-based Wells.
 */
class IWell : public ISourceSink {
public:
    int i, j, k_min, k_max;
    double q_total;
    bool is_injector;

    IWell(int i_v, int j_v, int kmin, int kmax, double q, bool inj)
        : i(i_v), j(j_v), k_min(kmin), k_max(kmax), q_total(q), is_injector(inj) {}
    
    virtual double get_q_water(const top::IState& state) const = 0;
};

/**
 * @brief ConstantRateWell for single-phase flow.
 */
class ConstantRateWell : public IWell {
public:
    using IndexFunc = std::function<int(int, int, int)>;
    IndexFunc get_idx;
    double scale_factor;

    ConstantRateWell(int i_v, int j_v, int kmin, int kmax, double q, double scale, IndexFunc idx)
        : IWell(i_v, j_v, kmin, kmax, q, false), scale_factor(scale), get_idx(idx) {}

    void apply(Vector& residual, Matrix* jacobian, const top::IState& state, double dt,
               std::vector<SparseMatrix::Entry>* sparse_entries = nullptr) override {
        int num_layers = k_max - k_min + 1;
        double q_scaled = (q_total * scale_factor) / num_layers;

        for (int k = k_min; k <= k_max; ++k) {
            int c = get_idx(i, j, k);
            residual[c] -= q_scaled;
        }
    }

    double get_q_water(const top::IState& state) const override { return 0.0; }
};

/**
 * @brief Base class for Dual-Phase Wells with Water/Oil physics.
 */
class DualPhaseWell : public IWell {
public:
    using RelPermFunc = std::function<void(double, double&, double&)>;
    using IndexFunc = std::function<int(int, int)>;
    using SwFunc = std::function<double(int, int)>;
    
    RelPermFunc get_rel_perm;
    IndexFunc get_idx;
    SwFunc get_sw; // Added to decouple from State class
    double mu_w, mu_o;

    DualPhaseWell(int i_v, int j_v, double q, bool inj, RelPermFunc rp, IndexFunc idx, SwFunc sw_f, double mw, double mo)
        : IWell(i_v, j_v, 0, 0, q, inj), get_rel_perm(rp), get_idx(idx), get_sw(sw_f), mu_w(mw), mu_o(mo) {}

    void apply(Vector& residual, Matrix* jacobian, const top::IState& state_raw, double dt, 
               std::vector<SparseMatrix::Entry>* sparse_entries = nullptr) override {}
    
    double get_q_water(const top::IState& state_raw) const override {
        if (is_injector) return q_total;
        
        double krw, kro;
        get_rel_perm(get_sw(i, j), krw, kro);
        double fw = (krw / mu_w) / ((krw / mu_w) + (kro / mu_o));
        return q_total * fw;
    }
};

/**
 * @brief Specialized well for 2D Dual Phase Reservoir (Fully Implicit support).
 */
class ReservoirWellDual2D : public DualPhaseWell {
public:
    ReservoirWellDual2D(int i_v, int j_v, double q, bool inj, RelPermFunc rp, IndexFunc idx, SwFunc sw_f, double mw, double mo)
        : DualPhaseWell(i_v, j_v, q, inj, rp, idx, sw_f, mw, mo) {}

    void apply(Vector& residual, Matrix* jacobian, const top::IState& state_raw, double dt,
               std::vector<SparseMatrix::Entry>* sparse_entries = nullptr) override {
        int c = get_idx(i, j);
        int r_w = 2 * c;
        int r_o = 2 * c + 1;

        if (is_injector) {
            residual[r_w] -= q_total; // Inject water
        } else {
            double krw, kro;
            double sw_val = get_sw(i, j);
            get_rel_perm(sw_val, krw, kro);
            double lam_w = krw / mu_w;
            double lam_o = kro / mu_o;
            double lam_t = lam_w + lam_o;
            double fw = lam_w / lam_t;
            
            residual[r_w] -= q_total * fw;
            residual[r_o] -= q_total * (1.0 - fw);

            auto apply_jac = [&](int r, int col, double val) {
                if (jacobian) (*jacobian)[r][col] += val;
                if (sparse_entries) sparse_entries->push_back({r, col, val});
            };

            if (jacobian || sparse_entries) {
                double sw_res = 0.2, so_res = 0.2; 
                double den_sw = 1.0 - sw_res - so_res;
                double swe = (sw_val - sw_res) / den_sw;
                double dkrw_dsw = (sw_val > sw_res && sw_val < (1.0 - so_res)) ? (2.0 * swe / den_sw) : 0.0;
                double dkro_dsw = (sw_val > sw_res && sw_val < (1.0 - so_res)) ? (-2.0 * (1.0 - swe) / den_sw) : 0.0;
                
                double dlamw_dsw = dkrw_dsw / mu_w;
                double dlamo_dsw = dkro_dsw / mu_o;
                double dfw_dsw = (dlamw_dsw * lam_o - lam_w * dlamo_dsw) / (lam_t * lam_t);
                
                apply_jac(r_w, r_w + 1, -q_total * dfw_dsw);
                apply_jac(r_o, r_w + 1, q_total * dfw_dsw);
            }
        }
    }
};

/**
 * @brief Specialized well for 3D Black Oil Reservoir (3-Phase support).
 */
class ReservoirWellBlackOil3D : public ISourceSink {
public:
    using RelPermFunc = std::function<void(double, double, double&, double&, double&)>;
    using IndexFunc = std::function<int(int, int, int)>;
    using VarFunc = std::function<double(int, int, int, int)>; // (i,j,k, var_idx)
    
    int i, j, k_min, k_max;
    double q_total;
    bool is_injector;
    
    RelPermFunc get_rel_perm;
    IndexFunc get_idx;
    VarFunc get_var; // Access P, Sw, Sg
    
    double mu_w, mu_o, mu_g;
    
    ReservoirWellBlackOil3D(int i_v, int j_v, int kmin, int kmax, double q, bool inj, 
                           RelPermFunc rp, IndexFunc idx, VarFunc var_f,
                           double mw, double mo, double mg)
        : i(i_v), j(j_v), k_min(kmin), k_max(kmax), q_total(q), is_injector(inj), 
          get_rel_perm(rp), get_idx(idx), get_var(var_f), 
          mu_w(mw), mu_o(mo), mu_g(mg) {}

    void apply(Vector& residual, Matrix* jacobian, const top::IState& state_raw, double dt,
               std::vector<SparseMatrix::Entry>* sparse_entries = nullptr) override {
        int num_layers = k_max - k_min + 1;
        double q_per_layer = q_total / num_layers;

        for (int k = k_min; k <= k_max; ++k) {
            int c = get_idx(i, j, k);
            int r_w = 3 * c;
            int r_o = 3 * c + 1;
            int r_g = 3 * c + 2;

            if (is_injector) {
                residual[r_w] -= q_per_layer; // Inject water
            } else {
                double krw, krog, krg;
                double sw_val = get_var(i, j, k, 1);
                double sg_val = get_var(i, j, k, 2);
                get_rel_perm(sw_val, sg_val, krw, krog, krg);
                
                double lam_w = krw / mu_w; 
                double lam_o = krog / mu_o;
                double lam_g = krg / mu_g;
                double lam_t = lam_w + lam_o + lam_g;
                if (lam_t < 1e-12) lam_t = 1e-12; // Safeguard
                
                residual[r_w] -= q_per_layer * (lam_w / lam_t);
                residual[r_o] -= q_per_layer * (lam_o / lam_t);
                residual[r_g] -= q_per_layer * (lam_g / lam_t);
            }
        }
    }
};

} // namespace mod
