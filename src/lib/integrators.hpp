#pragma once
#include "lib/base.hpp"
#include "lib/solvers.hpp"
#include <stdexcept>

namespace numerical_methods {

/**
 * @brief Forward Euler (Explicit) Integrator
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
 * @brief Implicit Euler Integrator
 * Requires a specialized LinearTridiagonalSolver for 1D systems.
 */
class ImplicitEulerIntegrator : public ITimeIntegrator {
public:
    void step(const IModel& model, IState& state, ISolver* solver, double dt) override {
        auto state_old = state.clone();

        if (auto* tri_solver = dynamic_cast<LinearTridiagonalSolver*>(solver)) {
            Vector R = model.build_residual(state, *state_old, dt);
            Matrix J = model.build_jacobian(state, dt);
            Vector delta = tri_solver->solve_system(J, R);
            state.update(delta);
        } 
        else if (auto* cg_solver = dynamic_cast<ConjugateGradientSolver*>(solver)) {
            // Setup matrix-vector operator for J * v
            auto apply_A = [&](const Vector& v) {
                return model.apply_jacobian(state, v, dt);
            };
            Vector R = model.build_residual(state, *state_old, dt);
            Vector guess(R.size(), 0.0);
            
            // Solve J * delta = -R
            Vector delta = cg_solver->solve_iterative(apply_A, scale(R, -1.0), guess);
            state.update(delta);
        }
        else {
            throw std::runtime_error("ImplicitEulerIntegrator requires a supported solver (Tridiagonal or CG).");
        }
    }

private:
    Vector scale(const Vector& v, double s) {
        Vector res = v;
        for (auto& val : res) val *= s;
        return res;
    }
};

/**
 * @brief Runge-Kutta 4 (Explicit) Integrator
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

} // namespace numerical_methods
