#pragma once
#include "lib/modules.hpp"
#include "lib/solvers.hpp"
#include "modules/reservoir/well.hpp"
#include <vector>
#include <algorithm>

namespace num {
using namespace top;
using namespace mod;

/**
 * @brief Reservoir IMPES Integrator (Implicit Pressure, Explicit Saturation).
 * 
 * Fits With:
 * - Two-phase or three-phase reservoir flow (e.g., Oil-Water, Black Oil).
 * - Decoupling pressure and saturation to reduce computational cost compared to Fully Implicit.
 * 
 * Best Case:
 * - Systems where pressure changes rapidly but saturation fronts move slowly.
 * - Pressure is solved implicitly (large dt stability), saturation is solved explicitly.
 * 
 * Worst Case:
 * - High flow rates or small grid cells where the explicit saturation update is limited by a strict CFL condition.
 * - Non-convergent pressure solutions if the Jacobian is not well-estimated.
 * 
 * Sample Input:
 * @tparam TState The specialized reservoir state class (e.g., Reservoir2DState).
 * @tparam TModel The specialized reservoir model class (e.g., DualPhaseModel).
 * @param solver Must be ConjugateGradientSolver for the implicit pressure stage.
 */
template<typename TState, typename TModel>
class ReservoirIMPESIntegrator : public ITimeIntegrator {
public:
    void step(const IModel& model_raw, IState& state_raw, ISolver* solver, double dt) override {
        auto& state = dynamic_cast<TState&>(state_raw);
        const auto& model = dynamic_cast<const TModel&>(model_raw);
        auto* cg_solver = dynamic_cast<ConjugateGradientSolver*>(solver);
        
        if (!cg_solver) throw std::runtime_error("IMPES requires ConjugateGradientSolver.");

        int nx = state.spatial.nx, ny = state.spatial.ny;
        double dx = state.spatial.dx, dy = state.spatial.dy;
        double unit_conv = 0.001127;
        
        // --- 1. SOLVE PRESSURE IMPLICITLY ---
        auto apply_A = [&](const Vector& v) {
            return model.apply_jacobian(state, v, dt);
        };
        
        Vector R_raw = model.build_residual(state, state, dt);
        Vector R_p(state.pressures.size(), 0.0);
        
        // Sum water and oil equations for pressure residual if coupled
        if (R_raw.size() == 2 * state.pressures.size()) {
            for (size_t i = 0; i < state.pressures.size(); ++i) {
                R_p[i] = R_raw[2 * i] + R_raw[2 * i + 1];
            }
        } else {
            R_p = R_raw;
        }

        Vector delta = cg_solver->solve_iterative(apply_A, scale(R_p, -1.0), Vector(R_p.size(), 0.0));
        state.update(delta);
        
        // --- 2. UPDATE SATURATION EXPLICITLY ---
        Vector sw_new = state.water_saturations;
        double pore_vol = model.phi * (dx * dy * model.h);

        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                int cur = state.idx(i, j);
                double water_flux_sum = 0.0;
                
                auto add_water_flux = [&](int ni, int nj, double dist) {
                    int neighbor = state.idx(ni, nj);
                    double p_cur = state.pressures[cur];
                    double p_neigh = state.pressures[neighbor];
                    
                    double sw_up = (p_cur > p_neigh) ? state.water_saturations[cur] : state.water_saturations[neighbor];
                    double krw, kro;
                    model.get_rel_perm(sw_up, krw, kro);
                    
                    double lambda_w = krw / model.mu_w;
                    double T_w = unit_conv * model.k_abs * lambda_w * (model.h * (dist == dx ? dy : dx)) / dist;
                    water_flux_sum += T_w * (p_neigh - p_cur);
                };

                if (i > 0) add_water_flux(i - 1, j, dx);
                if (i < nx - 1) add_water_flux(i + 1, j, dx);
                if (j > 0) add_water_flux(i, j - 1, dy);
                if (j < ny - 1) add_water_flux(i, j + 1, dy);

                double qw = 0.0;
                for (const auto& s : model.sources) {
                    auto w = std::dynamic_pointer_cast<mod::IWell>(s);
                    if (w && w->i == i && w->j == j) {
                        qw = w->get_q_water(state);
                    }
                }

                double delta_sw = (dt * (water_flux_sum + qw) * 5.615) / (24.0 * pore_vol);
                sw_new[cur] += delta_sw;
                sw_new[cur] = std::max(0.0, std::min(1.0, sw_new[cur]));
            }
        }
        state.water_saturations = sw_new;
    }

private:
    Vector scale(const Vector& v, double s) {
        Vector res = v;
        for (auto& val : res) val *= s;
        return res;
    }
};

} // namespace num
