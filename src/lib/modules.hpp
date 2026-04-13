/**
 * @file modules.hpp
 * @brief Architectural Contract between Physics & Numerics.
 * 
 * Objective:
 * This file defines the core structural abstractions (IState, IModel, ISolver, 
 * ITimeIntegrator) that enable the entire simulation suite to function.
 * 
 * Architectural Rationale:
 * Segregating these interfaces into a dedicated file provides a stable 
 * "contractual layer" that decouples WHAT we are simulating (Physics) 
 * from HOW we compute it (Numerics & Orchestration). 
 * 
 * Strategic Importance:
 * Without this segregation, physics modules would be tightly coupled to 
 * specific numeric solvers, preventing code reuse. This file enables 
 * "Plug-and-Play" extensibility, allowing a developer to add a new 
 * reservoir model without touching the time-stepping or linear algebra engines.
 */
#pragma once
#include <vector>
#include <memory>
#include <iostream>
#include <functional>

#include "lib/sparse.hpp"

namespace top {

// Placeholders for your actual math src 
using Vector = std::vector<double>; 
using Matrix = std::vector<std::vector<double>>; 
using SparseMatrix = num::SparseMatrix;

/**
 * @brief Abstract Base Class for Simulation State.
 * 
 * Each physics module (Heat, Reservoir, etc.) must implement its own State
 * to hold primary variables (Pressure, Temperature, Saturation).
 * 
 * Usage:
 * - Managed by ITimeIntegrator for temporal updates.
 * - Used by IModel to build residuals and Jacobians.
 */
class IState {
public:
    virtual ~IState() = default; // Crucial for polymorphic deletion
    
    /**
     * @brief Update state variables using a calculated delta vector.
     * @param delta The vector returned by a solver (typically x_{n+1} = x_n + delta).
     */
    virtual void update(const Vector& delta) = 0;
    
    /**
     * @brief Deep copy of the state object.
     * @return unique_ptr to the new clone.
     */
    virtual std::unique_ptr<IState> clone() const = 0; 
};

/**
 * @brief Abstract Base Class for Physical Models/Equations.
 * 
 * Defines the physics of the system through Residuals (R) and Jacobians (J).
 * 
 * Implementation Pattern:
 * - For Explicit methods: Implement evaluate_rhs().
 * - For Implicit methods: Implement build_residual() and build_jacobian().
 */
class IModel {
public:
    virtual ~IModel() = default;

    /**
     * @brief Evaluate f(u) for explicit schemes: du/dt = f(u).
     */
    virtual Vector evaluate_rhs(const IState& state) const = 0;
    
    /**
     * @brief Build the non-linear residual R = (u_new - u_old)/dt - f(u_new).
     */
    virtual Vector build_residual(const IState& state, const IState& state_old, double dt) const = 0;
    
    /**
     * @brief Build the Jacobian matrix dR/du. (Dense version, deprecated for large systems).
     */
    virtual Matrix build_jacobian(const IState& state, double dt) const = 0;

    /**
     * @brief Build the Sparse Jacobian matrix dR/du.
     */
    virtual SparseMatrix build_sparse_jacobian(const IState& state, double dt) const {
        return SparseMatrix(0, 0); // Default placeholder
    }

    /**
     * @brief Optional: Fast Matrix-Vector product J * v without full assembly.
     */
    virtual Vector apply_jacobian(const IState& state, const Vector& v, double dt) const {
        return {}; // Default no-op
    }
};

class ISolver {
public:
    virtual ~ISolver() = default;

    // Returns the delta vector to update the state
    virtual Vector solve(const IModel& model, const IState& state, double dt) = 0;
};

/**
 * @brief Base Class for Time-Stepping Strategies.
 */
class ITimeIntegrator {
public:
    virtual ~ITimeIntegrator() = default;

    /**
     * @brief Perform a single time step on the model/state.
     * @param model physics to integrate.
     * @param state state to update (modified in-place).
     * @param solver linear or non-linear solver to use if implicit.
     * @param dt time step size.
     */
    virtual void step(const IModel& model, IState& state, ISolver* solver, double dt) = 0;
};

/**
 * @brief High-level Orchestrator for physics simulations.
 * 
 * Uses Dependency Injection to combine Model, State, Solver, and Integrator.
 */
class StandardSimulator {
private:
    std::shared_ptr<IModel> model;
    std::shared_ptr<IState> state;
    std::shared_ptr<ISolver> solver;
    std::shared_ptr<ITimeIntegrator> integrator;

public:
    using StepCallback = std::function<void(double t, const IState& state)>;

    /**
     * @brief Constructor for the simulation environment.
     */
    StandardSimulator(std::shared_ptr<IModel> m, 
                      std::shared_ptr<IState> s, 
                      std::shared_ptr<ISolver> sol, 
                      std::shared_ptr<ITimeIntegrator> integ)
        : model(std::move(m)), state(std::move(s)), 
          solver(std::move(sol)), integrator(std::move(integ)) {}

    /**
     * @brief Run the simulation from t=0 to t_end.
     * @param on_step optional callback for telemetry or data export.
     */
    void run(double t_end, double dt, StepCallback on_step = nullptr) {
        double t = 0.0;
        int step_count = 0;
        
        if (on_step) on_step(t, *state);

        while (t < t_end) {
            integrator->step(*model, *state, solver.get(), dt);
            t += dt;
            step_count++;
            
            if (on_step) on_step(t, *state);
        }
    }
};

} // namespace top