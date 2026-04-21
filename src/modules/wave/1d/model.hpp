#pragma once
#include "lib/interfaces.hpp"
#include "state.hpp"
#include "lib/discretization.hpp"
 
namespace mod {
using namespace top;
namespace wave {

/**
 * @brief 1D Wave Physical Model (Properties).
 */
class Wave1DModel : public IModel {
public:
    std::shared_ptr<num::discretization::Conductance1D> cond;
    Vector storage_coeff;

    Wave1DModel(std::shared_ptr<num::discretization::Conductance1D> c, const Vector& storage) 
        : cond(c), storage_coeff(storage) {}

    double get_tolerance() const override { return 1e-4; }

    Vector get_accumulation_weights(const IGrid& grid, const IState& state) const override {
        size_t n = storage_coeff.size();
        Vector weights(2 * n);
        for (size_t i = 0; i < n; ++i) {
            weights[i] = 1.0;          // u accumulation weight
            weights[i + n] = storage_coeff[i]; // v accumulation weight (mass)
        }
        return weights;
    }
};

/**
 * @brief 1D Wave FVM Discretizer (Numerical Assembly).
 * Assembles the first-order system:
 *   du/dt = v
 *   m*dv/dt = Sum(Flux)
 */
class Wave1DDiscretizer : public IDiscretizer {
public:
    void assemble_jacobian(const IGrid& grid, const IModel& model, const IState& state, SparseMatrix& J) const override {
        const auto& w_model = static_cast<const Wave1DModel&>(model);
        const auto& w_state = static_cast<const Wave1DState&>(state);
        size_t n = w_state.u.size();

        if (J.rows != 2 * (int)n) J = SparseMatrix(2 * n, 2 * n);
        J.triplets.clear();

        // 1. du/dt = v part
        // J_uv = -I
        for (int i = 0; i < (int)n; ++i) {
            J.triplets.push_back({i, i + (int)n, -1.0});
        }

        // 2. m*dv/dt = Laplacian(u) part
        // J_vu = -Laplacian
        for (int i = 1; i < (int)n - 1; ++i) {
            double t_prev = w_model.cond->T[i-1];
            double t_next = w_model.cond->T[i];
            
            int cur_v = i + (int)n;
            J.triplets.push_back({cur_v, i - 1, -t_prev});
            J.triplets.push_back({cur_v, i, t_prev + t_next});
            J.triplets.push_back({cur_v, i + 1, -t_next});
        }
    }

    void assemble_residual(const IGrid& grid, const IModel& model, const IState& state, Vector& R) const override {
        const auto& w_model = static_cast<const Wave1DModel&>(model);
        const auto& w_state = static_cast<const Wave1DState&>(state);
        size_t n = w_state.u.size();

        // 1. R_u = -v
        #pragma omp parallel for
        for (int i = 0; i < (int)n; ++i) {
            R[i] = -w_state.v[i];
        }

        // 2. R_v = -NetFlux(u)
        #pragma omp parallel for
        for (int i = 1; i < (int)n - 1; ++i) {
            double net_flux = w_model.cond->T[i-1] * (w_state.u[i-1] - w_state.u[i]) +
                              w_model.cond->T[i]   * (w_state.u[i+1] - w_state.u[i]);
            R[i + n] = -net_flux;
        }
    }

    void apply_boundary_conditions(const IGrid& grid, const IModel& model, const IState& state, SparseMatrix& J, Vector& R) const override {
        const auto& w_state = static_cast<const Wave1DState&>(state);
        size_t n = w_state.u.size();

        // Dirichlet BCs for u (fixed displacement)
        R[0] = w_state.u[0]; // Assuming 0 for now or set from model
        R[n-1] = w_state.u[n-1];
        J.triplets.push_back({0, 0, 1.0});
        J.triplets.push_back({(int)n-1, (int)n-1, 1.0});

        // Fixed velocity at boundaries
        R[n] = w_state.v[0];
        R[2*n-1] = w_state.v[n-1];
        J.triplets.push_back({(int)n, (int)n, 1.0});
        J.triplets.push_back({(int)2*n-1, (int)2*n-1, 1.0});
    }
};

} // namespace wave
} // namespace mod
