/**
 * @file engine_infra.hpp
 * @brief Default implementations for new engine components.
 */
#pragma once
#include "interfaces.hpp"
#include <cmath>
#include <numeric>

namespace top {

/**
 * @brief Default Serial Implementation of ParallelManager.
 * Used for single-core execution.
 */
class SerialParallelManager : public IParallelManager {
public:
    void sync_ghost_cells(IState& state) const override {
        // No-op in serial
    }

    double get_global_norm(const Vector& r) const override {
        double local_norm = 0.0;
        for (double val : r) local_norm += val * val;
        return std::sqrt(local_norm);
    }
};

/**
 * @brief Minimal Serial Solver Wrapper (Stub).
 * Each concrete solver implementation will provide optimized solve logic.
 */
class SequentialSolver : public ISolver {
public:
    Vector solve(const SparseMatrix& A, const Vector& b) override {
        // This is a stub. Real solvers come from solvers.hpp
        return Vector(b.size(), 0.0);
    }
};

} // namespace top
