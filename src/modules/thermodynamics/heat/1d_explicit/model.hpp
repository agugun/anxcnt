#pragma once
#include "lib/interfaces.hpp"
#include "state.hpp"

namespace mod {

class Heat1DExplicitModel : public IModel {
private:
    std::shared_ptr<num::discretization::Conductance1D> cond;
    std::vector<double> storage_coeff;

public:
    Heat1DExplicitModel(std::shared_ptr<num::discretization::Conductance1D> c, const std::vector<double>& storage)
        : cond(c), storage_coeff(storage) {}

    void build_residual(const IState& st_new, const IState& st_old, double dt, std::vector<double>& R) const {
        const auto& s_new = static_cast<const Heat1DExplicitState&>(st_new);
        const auto& s_old = static_cast<const Heat1DExplicitState&>(st_old);
        size_t n = s_old.temperatures.size();

        #pragma omp parallel for
        for (int i = 0; i < (int)n; ++i) {
            double net_flux = 0.0;
            if (i > 0) net_flux += cond->T[i-1] * (s_old.temperatures[i-1] - s_old.temperatures[i]);
            if (i < (int)n - 1) net_flux += cond->T[i] * (s_old.temperatures[i+1] - s_old.temperatures[i]);

            R[i] = -(dt * net_flux) / storage_coeff[i];
        }
    }

    double get_tolerance() const override { return 1e-6; }

    Vector build_capacity(const IGrid& grd, const IState& st) const override {
        return storage_coeff;
    }
};

class Heat1DExplicitDiscretizer : public IDiscretizer {
public:
    void build_residual(const IGrid& grd, const IModel& mdl, const IState& st, Vector& R) const override {
        // For explicit, the mdl handles the RHS calculation.
    }

    void build_jacobian(const IGrid& grd, const IModel& mdl, const IState& st, SparseMatrix& J) const override {
        // Explicit Jacobian is not used in the same way.
    }

    void apply_bc(const IGrid& grd, const IModel& mdl, const IState& st, SparseMatrix& J, Vector& R) const override {
        // Boundary conditions are handled by the mdl.
    }
};

} // namespace mod
