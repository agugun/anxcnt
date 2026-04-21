/**
 * @file interfaces.hpp
 * @brief Centralized Component-Based Architecture for Numerical Simulations.
 */
#pragma once
#include <vector>
#include <memory>
#include <string>
#include <functional>
#include "lib/sparse.hpp"

namespace top {

// Global types for consistency
using Vector = std::vector<double>;
using Matrix = std::vector<std::vector<double>>;
using SparseMatrix = num::SparseMatrix;

/**
 * @brief Handles topology and geometry.
 * Who is next to whom, and how big are they?
 */
class IGrid {
public:
    virtual ~IGrid() = default;
    virtual size_t get_total_cells() const = 0;
};

/**
 * @brief Handles variable memory and physical bounds.
 */
class IState {
public:
    virtual ~IState() = default;
    virtual void apply_update(const std::vector<double>& delta) = 0;
    virtual std::vector<double> to_vector() const = 0;
    virtual std::unique_ptr<IState> clone() const = 0;
};

/**
 * @brief Generic property lookups (PVT, Rock, Thermal, etc).
 */
class IPropertyManager {
public:
    virtual ~IPropertyManager() = default;
};

/**
 * @brief Translates continuous physics into discrete terms.
 */
class IModel {
public:
    virtual ~IModel() = default;
    virtual double get_tolerance() const = 0;
    
    /**
     * @brief Returns the diagonal weights for the accumulation term (e.g., C in C*du/dt).
     */
    virtual Vector get_accumulation_weights(const IGrid& grid, const IState& state) const = 0;
};

/**
 * @brief Numerical Mapping: Translates Grid + Model into Sparse Matrices.
 * Covers FDM, FVM, and FEM assembly logic.
 */
class IDiscretizer {
public:
    virtual ~IDiscretizer() = default;
    virtual void assemble_jacobian(const IGrid& grid, const IModel& model, const IState& state, SparseMatrix& J) const = 0;
    virtual void assemble_residual(const IGrid& grid, const IModel& model, const IState& state, Vector& R) const = 0;
    virtual void apply_boundary_conditions(const IGrid& grid, const IModel& model, const IState& state, SparseMatrix& J, Vector& R) const = 0;
};

/**
 * @brief External Source/Sink or Inner Boundary Terms.
 */
class ISourceSink {
public:
    virtual ~ISourceSink() = default;
    virtual void assemble_terms(const IState& state, SparseMatrix& J, Vector& R) const = 0;
};

/**
 * @brief Temporal Strategy: Manages accumulation terms (du/dt).
 */
class ITimeIntegrator {
public:
    virtual ~ITimeIntegrator() = default;
    virtual double get_next_timestep(const IState& state, double t) const = 0;
    virtual void add_accumulation(const IGrid& grid, const IModel& model, SparseMatrix& J, Vector& R, const IState& state_new, const IState& state_old, double dt) const = 0;
};

/**
 * @brief Pure Algebraic Solver (Ax = b).
 */
class ISolver {
public:
    virtual ~ISolver() = default;
    virtual Vector solve(const SparseMatrix& A, const Vector& b) = 0;
};

/**
 * @brief Parallel Orchestrator (MPI/Halo-sync and Global Reductions).
 */
class IParallelManager {
public:
    virtual ~IParallelManager() = default;
    virtual void sync_ghost_cells(IState& state) const = 0;
    virtual double get_global_norm(const Vector& r) const = 0;
};

/**
 * @brief Non-linear Orchestrator (Orchestrates assembly and solve steps).
 */
class ILinearizer {
public:
    virtual ~ILinearizer() = default;
    virtual std::unique_ptr<IState> solve_timestep(
        const IState& state_n, double dt, 
        const IGrid& grid, const IModel& model, const IDiscretizer& discretizer,
        const ITimeIntegrator& timer, ISolver& solver, const IParallelManager& pm,
        const std::vector<std::shared_ptr<ISourceSink>>& sources) = 0;
};

/**
 * @brief Telemetry and Output (Handles logging, VTK/CSV exporting, etc).
 */
class IObserver {
public:
    virtual ~IObserver() = default;
    virtual void on_simulation_start(const IGrid& grid) {}
    virtual void on_step_complete(double t, int step, const IState& state) = 0;
    virtual void on_simulation_end() {}
};

} // namespace top
