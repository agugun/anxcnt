#pragma once
#include "lib/interfaces.hpp"
#include "state.hpp"
#include "model.hpp"

namespace mod {
using namespace top;
namespace pressure {

/**
 * @brief 1D Pressure FVM Discretizer (Numerical Assembly).
 */
class Pressure1DDiscretizer : public IDiscretizer {
public:
    void assemble_jacobian(const IGrid& grid, const IModel& model, const IState& state, SparseMatrix& J) const override {
        const auto& p_model = static_cast<const Pressure1DModel&>(model);
        const auto& p_state = static_cast<const Pressure1DState&>(state);
        size_t n = p_state.pressures.size();

        if (J.rows != n) J = SparseMatrix(n, n);
        J.triplets.clear();

        for (int i = 1; i < (int)n - 1; ++i) {
            double t_prev = p_model.cond->T[i-1];
            double t_next = p_model.cond->T[i];
            
            J.triplets.push_back({i, i-1, -t_prev});
            J.triplets.push_back({i, i, t_prev + t_next});
            J.triplets.push_back({i, i+1, -t_next});
        }
    }

    void assemble_residual(const IGrid& grid, const IModel& model, const IState& state, Vector& R) const override {
        const auto& p_model = static_cast<const Pressure1DModel&>(model);
        const auto& p_state = static_cast<const Pressure1DState&>(state);
        size_t n = p_state.pressures.size();

        #pragma omp parallel for
        for (int i = 1; i < (int)n - 1; ++i) {
            double net_flux = p_model.cond->T[i-1] * (p_state.pressures[i-1] - p_state.pressures[i]) +
                              p_model.cond->T[i]   * (p_state.pressures[i+1] - p_state.pressures[i]);
            R[i] = -net_flux;
        }
    }

    void apply_boundary_conditions(const IGrid& grid, const IModel& model, const IState& state, SparseMatrix& J, Vector& R) const override {
        const auto& p_model = static_cast<const Pressure1DModel&>(model);
        const auto& p_state = static_cast<const Pressure1DState&>(state);
        size_t n = R.size();

        // Dirichlet BCs
        R[0] = p_state.pressures[0] - p_model.p_left;
        R[n-1] = p_state.pressures[n-1] - p_model.p_right;

        J.triplets.push_back({0, 0, 1.0});
        J.triplets.push_back({(int)n-1, (int)n-1, 1.0});
    }
};

} // namespace pressure
} // namespace mod
