#pragma once
#include "lib/modules.hpp"
#include "state.hpp"
#include "lib/operators.hpp"
#include "modules/reservoir/well.hpp"
#include <cmath>
#include <vector>

namespace mod {
using namespace top;
namespace reservoir {

/**
 * @brief Classical 1D Pressure Diffusivity Reservoir Model.
 * 
 * Governing Equation:
 *   d^2P/dx^2 = (1 / eta) * dP/dt
 * 
 * Physical Logic:
 *   Models single-phase fluid flow in a porous medium. This is the 1D line-source 
 *   solution equivalent.
 *   - eta: Hydraulic diffusivity = 0.0002637 * k / (phi * mu * ct) [ft^2/hr]
 * 
 * Discretization:
 *   - Spatial: 2nd Order Central Difference for Laplacian.
 *   - Temporal: Backward Euler (Implicit) method.
 * 
 * Boundary Conditions:
 *   Default behavior in the Jacobian construction implements No-flow (Neumann) 
 *   boundaries (dP/dx = 0) at both ends by assuming p_{-1} = p_0 and p_N = p_{N-1}.
 */
class Reservoir1DModel : public IModel {
private:
    // Rock/Fluid Properties
    double k;      // permeability [mD]
    double phi;    // porosity [fraction]
    double mu;     // viscosity [cP]
    double ct;     // total compressibility [psi^-1]
    double B;      // formation volume factor [rb/stb]
    double area;   // cross-sectional area [ft^2]
    
    std::vector<std::shared_ptr<ISourceSink>> sources;

    double eta; // diffusivity [ft^2 / hr]

public:
    Reservoir1DModel(double k_val, double phi_val, double mu_val, double ct_val, 
                     double B_val, double area_val, const std::vector<std::shared_ptr<ISourceSink>>& sources_val)
        : k(k_val), phi(phi_val), mu(mu_val), ct(ct_val), 
          B(B_val), area(area_val), sources(sources_val) {
        
        // Oil field units constant: 0.0002637 for ft^2/hr
        eta = 0.0002637 * k / (phi * mu * ct);
    }

    Vector evaluate_rhs(const IState& state) const override {
        const auto& r_state = dynamic_cast<const Reservoir1DState&>(state);
        Vector rhs = mop::laplace_1d(r_state.pressures, r_state.spatial.dx);
        for (auto& v : rhs) v *= eta;
        
        for (auto& s : sources) {
            s->apply(rhs, nullptr, r_state, 0.0);
        }

        return rhs;
    }

    Vector build_residual(const IState& s_new, const IState& s_old, double dt) const override {
        const auto& r_new = dynamic_cast<const Reservoir1DState&>(s_new);
        const auto& r_old = dynamic_cast<const Reservoir1DState&>(s_old);
        
        Vector lap = mop::laplace_1d(r_new.pressures, r_new.spatial.dx);
        Vector res(r_new.pressures.size());
        for (size_t i = 0; i < res.size(); ++i) {
            res[i] = r_new.pressures[i] - r_old.pressures[i] - dt * eta * lap[i];
        }

        Vector well_contributions(res.size(), 0.0);
        for (auto& s : sources) {
            s->apply(well_contributions, nullptr, r_new, dt);
        }
        for (size_t i = 0; i < res.size(); ++i) {
            res[i] += dt * well_contributions[i];
        }

        return res;
    }

    Matrix build_jacobian(const IState& state, double dt) const override {
        const auto& r_state = dynamic_cast<const Reservoir1DState&>(state);
        size_t n = r_state.pressures.size();
        double dx = r_state.spatial.dx;
        double trans = eta * dt / (dx * dx);

        Matrix jac(n, Vector(n, 0.0));
        for (size_t i = 0; i < n; ++i) {
            jac[i][i] = 1.0 + 2.0 * trans; // Base diagonal
            
            if (i > 0) {
                jac[i][i-1] = -trans;
            } else {
                // Left No-flow: p_{-1} = p_0 => (p_1 - 2p_0 + p_0) terms
                jac[i][i] -= trans; 
            }
            
            if (i < n - 1) {
                jac[i][i+1] = -trans;
            } else {
                // Right No-flow: p_{n} = p_{n-1}
                jac[i][i] -= trans;
            }
        }
        return jac;
    }
};

} // namespace reservoir
} // namespace mod
