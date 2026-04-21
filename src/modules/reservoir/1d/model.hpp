#pragma once
#include "lib/interfaces.hpp"
#include "state.hpp"
#include "lib/discretization.hpp"
#include "modules/reservoir/well.hpp"

namespace mod {
using namespace top;
namespace reservoir {

/**
 * @brief 1D Reservoir Physical Model (Properties).
 */
class Reservoir1DModel : public IModel {
public:
    std::shared_ptr<num::discretization::Conductance1D> cond;
    Vector storage_coeff;
    std::vector<std::shared_ptr<ISourceSink>> wells;

    Reservoir1DModel(std::shared_ptr<num::discretization::Conductance1D> c, const Vector& storage, 
                     const std::vector<std::shared_ptr<ISourceSink>>& wells_val)
        : cond(c), storage_coeff(storage), wells(wells_val) {}

    double get_tolerance() const override { return 1e-4; }

    Vector get_accumulation_weights(const IGrid& grid, const IState& state) const override {
        return storage_coeff;
    }
    
    const std::vector<std::shared_ptr<ISourceSink>>& get_sources() const {
        return wells;
    }
};

/**
 * @brief 1D Reservoir FVM Discretizer.
 */
class Reservoir1DDiscretizer : public IDiscretizer {
public:
    void assemble_jacobian(const IGrid& grid, const IModel& model, const IState& state, SparseMatrix& J) const override {
        const auto& r_model = static_cast<const Reservoir1DModel&>(model);
        const auto& r_state = static_cast<const Reservoir1DState&>(state);
        size_t n = r_state.pressures.size();

        if (J.rows != (int)n) J = SparseMatrix(n, n);
        J.triplets.clear();

        for (int i = 0; i < (int)n; ++i) {
            double diag = 0.0;
            if (i > 0) {
                double t_prev = r_model.cond->T[i-1];
                diag += t_prev;
                J.triplets.push_back({i, i - 1, -t_prev});
            }
            if (i < (int)n - 1) {
                double t_next = r_model.cond->T[i];
                diag += t_next;
                J.triplets.push_back({i, i + 1, -t_next});
            }
            J.triplets.push_back({i, i, diag});
        }
    }

    void assemble_residual(const IGrid& grid, const IModel& model, const IState& state, Vector& R) const override {
        const auto& r_model = static_cast<const Reservoir1DModel&>(model);
        const auto& r_state = static_cast<const Reservoir1DState&>(state);
        size_t n = r_state.pressures.size();

        #pragma omp parallel for
        for (int i = 0; i < (int)n; ++i) {
            double net_flux = 0.0;
            if (i > 0)          net_flux += r_model.cond->T[i-1] * (r_state.pressures[i-1] - r_state.pressures[i]);
            if (i < (int)n - 1) net_flux += r_model.cond->T[i]   * (r_state.pressures[i+1] - r_state.pressures[i]);
            R[i] = -net_flux;
        }
    }

    void apply_boundary_conditions(const IGrid& grid, const IModel& model, const IState& state, SparseMatrix& J, Vector& R) const override {
        // Reservoir 1D usually assumes no-flow by default in flux assembly.
        // Dirichlet can be added here if needed.
    }
};

} // namespace reservoir
} // namespace mod
