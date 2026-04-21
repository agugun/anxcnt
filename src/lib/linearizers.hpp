/**
 * @file linearizers.hpp
 * @brief Strategies for hitting mathematical convergence within a time step.
 */
#pragma once
#include "interfaces.hpp"
#include <iostream>
#include <cmath>

namespace top {

/**
 * @brief Standard Newton-Raphson non-linear solver.
 */
class NewtonRaphson : public ILinearizer {
private:
    int max_iter;
    double tolerance;
    bool verbose;

public:
    NewtonRaphson(double tol = 1e-4, int max_it = 12, bool verb = true)
        : max_iter(max_it), tolerance(tol), verbose(verb) {}

    std::unique_ptr<IState> solve_timestep(
        const IState& state_n, double dt, 
        const IGrid& grid, const IModel& model, const IDiscretizer& discretizer,
        const ITimeIntegrator& timer, ISolver& solver, const IParallelManager& pm,
        const std::vector<std::shared_ptr<ISourceSink>>& sources) override 
    {
        auto state_guess = state_n.clone();
        bool converged = false;
        size_t n_cells = grid.get_total_cells();

        // SparseMatrix initialization logic
        SparseMatrix J(n_cells, n_cells);
        Vector R(n_cells, 0.0);

        for (int iter = 0; iter < max_iter; ++iter) {
            // 1. Parallel Sync (Update ghost/halo cells)
            pm.sync_ghost_cells(*state_guess);

            // Reset Residual and Jacobian triplets (Assembly logic)
            R.assign(n_cells, 0.0);
            J.clear();
            
            // 2. Spatial Assembly (Fluxes and properties)
            discretizer.assemble_jacobian(grid, model, *state_guess, J);
            discretizer.assemble_residual(grid, model, *state_guess, R);

            // 3. Source/Sink Assembly (Wells/Boundary Terms)
            for (auto& source : sources) {
                source->assemble_terms(*state_guess, J, R);
            }

            // 4. External Boundary Conditions
            discretizer.apply_boundary_conditions(grid, model, *state_guess, J, R);

            // 5. Temporal Accumulation (du/dt)
            timer.add_accumulation(grid, model, J, R, *state_guess, state_n, dt);

            // Finalize Matrix for Solve
            J.compress();

            // 6. Global Convergence Check
            double global_error = pm.get_global_norm(R);
            if (verbose) {
                std::cout << "    Newton Iter " << iter << ": Res Norm = " << global_error << std::endl;
            }

            if (global_error < tolerance) {
                converged = true;
                break;
            }

            // 7. Algebraic Solve: J * delta = -R
            // We pass -R to the solver
            Vector neg_R = R;
            for (auto& val : neg_R) val *= -1.0;
            
            Vector delta_x = solver.solve(J, neg_R);

            // 8. Safely Update State
            state_guess->apply_update(delta_x);
        }

        if (!converged && verbose) {
            std::cerr << "    [NewtonRaphson] Warning: Did not converge within " << max_iter << " iterations." << std::endl;
        }

        return state_guess;
    }
};

/**
 * @brief Explicit Linearizer for non-iterative steps (Explicit schemes).
 */
class ExplicitLinearizer : public ILinearizer {
public:
    std::unique_ptr<IState> solve_timestep(
        const IState& state_n, double dt, 
        const IGrid& grid, const IModel& model, const IDiscretizer& discretizer,
        const ITimeIntegrator& timer, ISolver& solver, const IParallelManager& pm,
        const std::vector<std::shared_ptr<ISourceSink>>& sources) override 
    {
        auto state_guess = state_n.clone();
        size_t n_cells = grid.get_total_cells();

        Vector R(n_cells, 0.0);
        
        // 1. Parallel Sync
        pm.sync_ghost_cells(*state_guess);

        // 2. Spatial Assembly (Fluxes only, no Jacobian needed for Pure Explicit)
        discretizer.assemble_residual(grid, model, state_n, R);

        // 3. Source/Sink Assembly 
        SparseMatrix J_dummy(0,0); // Dummy for explicit
        for (auto& source : sources) {
            source->assemble_terms(state_n, J_dummy, R);
        }

        // 4. Update state directly: state_n+1 = state_n + delta (where R already includes dt and storage effects)
        state_guess->apply_update(R);

        return state_guess;
    }
};

} // namespace top
