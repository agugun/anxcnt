#pragma once
#include "i_simulator.hpp"
#include "../model/state/field.hpp"
#include "../model/system_matrix.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace simulation {

/**
 * @brief Time evolution engine using Newton-Raphson method for nonlinear problems.
 *        Compatible with interactive steps (UI) and batch runs.
 */
class NewtonSimulator : public ISimulator {
public:
    NewtonSimulator() : max_iterations(10), tolerance(1e-6) {}

    void set_solver(std::shared_ptr<numerical_methods::solver::ISolver> s) override {
        solver = s;
    }

    void set_max_iterations(int max_iter) { max_iterations = max_iter; }
    void set_tolerance(double tol) { tolerance = tol; }

    /**
     * @brief Run a full batch time loop.
     */
    void run(model::IModel& model, model::state::IState& state, const std::map<std::string, double>& time_params) override {
        double dt = time_params.at("dt");
        double steps = time_params.at("steps");
        double current_time = 0.0;

        for (int i = 0; i < static_cast<int>(steps); ++i) {
            step(dt, model, state);
            current_time += dt;
        }
    }

    /**
     * @brief Execute one completely implicit Newton-Raphson timestep.
     */
    void step(double dt, model::IModel& model, model::state::IState& state) {
        if (!solver) throw std::runtime_error("Solver not set in NewtonSimulator.");

        // Downcast to access copy easily, or we can assume state provides copy mechanism.
        // For physics model we use Field. Let's just create raw fields.
        auto* p_curr = dynamic_cast<model::state::Field*>(&state);
        if(!p_curr) throw std::runtime_error("State is not a Field");
        
        model::state::Field p_prev = *p_curr;
        model.set_previous_state(p_prev);

        for (int iter = 0; iter < max_iterations; ++iter) {
            model::SystemMatrix sys;
            model.build_system(dt, *p_curr, sys);

            // Calculate residual norm (stored in sys.d)
            double residual_norm = 0.0;
            for (double val : sys.d) {
                residual_norm += val * val;
            }
            residual_norm = std::sqrt(residual_norm);

            if (residual_norm < tolerance) {
                // Converged
                return;
            }

            // Solve J * delta = r -> Linear solver handles this
            model::SystemMatrix J = sys; 
            // In our LinearSolver implementation it will directly replace internal points of p_curr.
            // But we want it to subtract (p^{k+1} = p^{k} - delta_p).
            // Let's create a temporary delta field.
            model::state::Field delta(p_curr->nx, p_curr->dx);
            solver->solve(J, delta);

            for (int i = 0; i < p_curr->nx - 2; ++i) {
                p_curr->values[i + 1] -= delta.values[i + 1];
            }
        }
        
        throw std::runtime_error("Newton's method failed to converge within the maximum number of iterations.");
    }

    void write_output(std::ostream& os, const model::state::IState& state, double current_time) override {
        const auto* field = dynamic_cast<const model::state::Field*>(&state);
        if (!field) return;

        // Tidy data format: Time, Index, X, Value
        const auto& data = field->get_data();
        for (int i = 0; i < field->nx; ++i) {
            os << current_time << "," 
               << i << "," 
               << field->get_coordinate(i) << "," 
               << data[i] << "\n";
        }
    }

private:
    std::shared_ptr<numerical_methods::solver::ISolver> solver;
    int max_iterations;
    double tolerance;
};

} // namespace simulation
