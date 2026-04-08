#pragma once
#include "i_solver.hpp"
#include "../../model/system_matrix.hpp"
#include "tridiagonal_solver.hpp"
#include "../../model/state/field.hpp"
#include <stdexcept>

namespace numerical_methods::solver {

class LinearSolver : public ISolver {
public:
    void solve(const model::SystemMatrix& sys, model::state::IState& state) override {
        // Downcast IState to Field
        auto* field = dynamic_cast<model::state::Field*>(&state);
        if (!field) {
            throw std::runtime_error("LinearSolver requires a Field state.");
        }

        // Solve tridiagonal system
        std::vector<double> sol = TridiagonalSolver::solve(sys.a, sys.b, sys.c, sys.d);

        // Update interior points
        for (size_t i = 0; i < sol.size(); ++i) {
            field->values[i + 1] = sol[i];
        }
    }
};

} // namespace numerical_methods::solver
