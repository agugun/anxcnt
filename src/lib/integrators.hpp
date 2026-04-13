/**
 * @file integrators.hpp
 * @brief Temporal Orchestration Engine.
 * 
 * Objective:
 * This file handles the logic for advancing simulation states through time
 * by coordinating between physical models and algebraic solvers.
 * 
 * Architectural Rationale:
 * Segregating integrators separates the "Flow of Time" from the "Calculus 
 * of Space." By isolating time-stepping patterns (e.g., Forward Euler, 
 * RK4, Implicit Euler) from the specific physics equations, we ensure 
 * that the system's temporal evolution remains independent of its 
 * physical constraints.
 * 
 * Strategic Importance:
 * It serves as the "Execution Glue" that binds Models and Solvers only 
 * during the runtime loop. This segregation allows a user to change 
 * the simulation accuracy (e.g., upgrading from 1st-order Euler to 
 * 4th-order Runge-Kutta) without changing the underlying PDE implementation.
 */
#pragma once
#include "lib/modules.hpp"
#include "lib/solvers.hpp"
#include <stdexcept>
#include <vector>

namespace num {
using namespace top;


/**
 * @brief Forward Euler (Explicit) Integrator.
 * 
 * Fits With:
 * - Simple linear or weakly non-linear ODEs where stability is not the primary constraint.
 * - Systems where dt satisfies the CFL (Courant-Friedrichs-Lewy) condition.
 * 
 * Best Case:
 * - Low-cost step evaluation. Accuracy is O(dt).
 * 
 * Worst Case:
 * - Stiff systems (requires extremely small dt to avoid divergence).
 * 
 * Sample Input:
 * @param model IModel implementation (requires evaluate_rhs).
 * @param state IState implementation (updated in-place).
 * @param solver Not used for explicit methods.
 * @param dt Time step size.
 */
class ForwardEulerIntegrator : public ITimeIntegrator {
public:
    void step(const IModel& model, IState& state, ISolver* solver, double dt) override {
        Vector rhs = model.evaluate_rhs(state);
        Vector delta(rhs.size());
        for (size_t i = 0; i < rhs.size(); ++i) {
            delta[i] = dt * rhs[i];
        }
        state.update(delta);
    }
};

/**
 * @brief Implicit Euler Integrator.
 * 
 * Fits With:
 * - Stiff systems where Forward Euler requires too small a time step.
 * - Diffusion problems (e.g., Heat Equation, Pressure Diffusivity).
 * 
 * Best Case:
 * - Unconditionally stable for linear problems. Allows large time steps.
 * 
 * Worst Case:
 * - First-order accurate O(dt). May cause numerical diffusion.
 * - Requires linear solver (CG or Tridiagonal) each step.
 */
class ImplicitEulerIntegrator : public ITimeIntegrator {
public:
    void step(const IModel& model, IState& state, ISolver* solver, double dt) override {
        if (!solver) {
            throw std::runtime_error("ImplicitEulerIntegrator requires a valid ISolver implementation.");
        }
        
        // Pass orchestration to the solver's unified interface
        Vector delta = solver->solve(model, state, dt);
        state.update(delta);
    }
};

/**
 * @brief Runge-Kutta 4 (Explicit) Integrator.
 * 
 * Fits With:
 * - Dynamic systems where high temporal accuracy is required (e.g., Wave Equation).
 * 
 * Best Case:
 * - Fourth-order accurate O(dt^4). Much more accurate than Euler for the same dt.
 * 
 * Worst Case:
 * - 4x more expensive per step than Forward Euler.
 * - Still limited by CFL stability constraints.
 */
class RungeKutta4Integrator : public ITimeIntegrator {
public:
    void step(const IModel& model, IState& state, ISolver* solver, double dt) override {
        // k1
        Vector k1 = model.evaluate_rhs(state);

        // k2
        auto s2 = state.clone();
        s2->update(scale(k1, dt * 0.5));
        Vector k2 = model.evaluate_rhs(*s2);

        // k3
        auto s3 = state.clone();
        s3->update(scale(k2, dt * 0.5));
        Vector k3 = model.evaluate_rhs(*s3);

        // k4
        auto s4 = state.clone();
        s4->update(scale(k3, dt));
        Vector k4 = model.evaluate_rhs(*s4);

        // Combine
        Vector delta(k1.size());
        for (size_t i = 0; i < k1.size(); ++i) {
            delta[i] = (dt / 6.0) * (k1[i] + 2.0 * (k2[i] + k3[i]) + k4[i]);
        }
        state.update(delta);
    }

private:
    Vector scale(const Vector& v, double s) {
        Vector res = v;
        for (auto& val : res) val *= s;
        return res;
    }
};

/**
 * @brief Fully Implicit (Non-linear) Integrator.
 * 
 * Fits With:
 * - Strongly non-linear coupled systems (e.g., Black Oil reservoir simulation).
 * 
 * Best Case:
 * - Extremely stable for highly non-linear residuals.
 * 
 * Worst Case:
 * - Most expensive. Requires multiple Newton iterations per time step.
 */
class FullyImplicitIntegrator : public ITimeIntegrator {
public:
    void step(const IModel& model, IState& state, ISolver* solver, double dt) override {
        if (!solver) {
            throw std::runtime_error("FullyImplicitIntegrator requires a valid ISolver implementation.");
        }
        
        Vector delta = solver->solve(model, state, dt);
        state.update(delta);
    }
};

} // namespace num
