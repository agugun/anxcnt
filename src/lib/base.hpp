#pragma once
#include <vector>
#include <memory>
#include <iostream>
#include <functional>

// Placeholders for your actual math src 
using Vector = std::vector<double>; 
using Matrix = std::vector<std::vector<double>>; 

class IState {
public:
    virtual ~IState() = default; // Crucial for polymorphic deletion
    
    // Updates the internal fields (e.g., pressure, saturation)
    virtual void update(const Vector& delta) = 0;
    
    // The Virtual Constructor idiom (replaces Python's copy)
    virtual std::unique_ptr<IState> clone() const = 0; 
};

class IModel {
public:
    virtual ~IModel() = default;

    // Const guarantees the model equations aren't altered during evaluation
    virtual Vector evaluate_rhs(const IState& state) const = 0;
    
    virtual Vector build_residual(const IState& state, const IState& state_old, double dt) const = 0;
    
    virtual Matrix build_jacobian(const IState& state, double dt) const = 0;

    // Optional: Matrix-vector product for iterative solvers (J * v)
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

class ITimeIntegrator {
public:
    virtual ~ITimeIntegrator() = default;

    // State is passed by reference so it can be updated in place
    virtual void step(const IModel& model, IState& state, ISolver* solver, double dt) = 0;
};

class StandardSimulator {
private:
    std::shared_ptr<IModel> model;
    std::shared_ptr<IState> state;
    std::shared_ptr<ISolver> solver;
    std::shared_ptr<ITimeIntegrator> integrator;

public:
    using StepCallback = std::function<void(double t, const IState& state)>;

    // Inject the dependencies via the constructor
    StandardSimulator(std::shared_ptr<IModel> m, 
                      std::shared_ptr<IState> s, 
                      std::shared_ptr<ISolver> sol, 
                      std::shared_ptr<ITimeIntegrator> integ)
        : model(std::move(m)), state(std::move(s)), 
          solver(std::move(sol)), integrator(std::move(integ)) {}

    void run(double t_end, double dt, StepCallback on_step = nullptr) {
        double t = 0.0;
        int step_count = 0;
        
        // Initial call if needed
        if (on_step) on_step(t, *state);

        while (t < t_end) {
            // The simulator orchestrates the step without knowing the math
            integrator->step(*model, *state, solver.get(), dt);
            t += dt;
            step_count++;
            
            if (on_step) on_step(t, *state);
        }
    }
};