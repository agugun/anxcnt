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

namespace geo {

/**
 * @brief Handles topology and geometry.
 */
class IGrid {
public:
    virtual ~IGrid() = default;
    virtual size_t get_total_cells() const = 0;
};

} // namespace geo

namespace mod {

// Global types for consistency
using Vector = std::vector<double>;
using Matrix = std::vector<std::vector<double>>;
using SparseMatrix = num::SparseMatrix;

// Backwards-compatible alias while geometry contracts move to namespace geo.
using IGrid = geo::IGrid;

/**
 * @brief Handles variable memory and physical bounds.
 */
class IState {
public:
    virtual ~IState() = default;
    virtual void update(const std::vector<double>& delta) = 0;
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
    virtual Vector build_capacity(const geo::IGrid& grd, const IState& st) const = 0;
};

/**
 * @brief Numerical Mapping: Translates Grid + Model into Sparse Matrices.
 */
class IDiscretizer {
public:
    virtual ~IDiscretizer() = default;
    virtual void build_jacobian(const geo::IGrid& grd, const IModel& mdl, const IState& st, SparseMatrix& J) const = 0;
    virtual void build_residual(const geo::IGrid& grd, const IModel& mdl, const IState& st, Vector& R) const = 0;
    virtual void apply_bc(const geo::IGrid& grd, const IModel& mdl, const IState& st, SparseMatrix& J, Vector& R) const = 0;
};

/**
 * @brief External Source/Sink or Inner Boundary Terms.
 */
class ISourceSink {
public:
    virtual ~ISourceSink() = default;
    virtual void assemble_terms(const IState& st, SparseMatrix& J, Vector& R) const = 0;
};

} // namespace mod

namespace utl {

/**
 * @brief Parallel orchestration and global reductions.
 */
class IParallelManager {
public:
    virtual ~IParallelManager() = default;
    virtual void sync_ghost_cells(mod::IState& st) const = 0;
    virtual double get_global_norm(const mod::Vector& r) const = 0;
};

/**
 * @brief Telemetry and Output (Handles logging, VTK/CSV exporting, etc).
 */
class IObserver {
public:
    virtual ~IObserver() = default;
    virtual void on_simulation_start(const geo::IGrid& grd) {}
    virtual void on_step_complete(double t, int step, const mod::IState& st) = 0;
    virtual void on_simulation_end() {}
};

} // namespace utl

namespace num {

/**
 * @brief Temporal Strategy: Manages accumulation terms (du/dt).
 */
class ITimeIntegrator {
public:
    virtual ~ITimeIntegrator() = default;
    virtual double compute_dt(const mod::IState& st, double t) const = 0;
    virtual void apply_temporal(const geo::IGrid& grd, const mod::IModel& mdl, mod::SparseMatrix& J, mod::Vector& R, const mod::IState& st_new, const mod::IState& st_old, double dt) const = 0;
};

/**
 * @brief Pure Algebraic Solver (Ax = b).
 */
class ISolver {
public:
    virtual ~ISolver() = default;
    virtual mod::Vector solve(const mod::SparseMatrix& A, const mod::Vector& b) = 0;
};

/**
 * @brief Non-linear Orchestrator (Orchestrates assembly and solve steps).
 */
class ILinearizer {
public:
    virtual ~ILinearizer() = default;

    virtual void set_sources(const std::vector<std::shared_ptr<mod::ISourceSink>>& sources) = 0;

    virtual std::unique_ptr<mod::IState> resolve(
        const mod::IState& st_n, double dt,
        const geo::IGrid& grd, const mod::IModel& mdl, const mod::IDiscretizer& discretizer,
        const ITimeIntegrator& timer, ISolver& solver, const utl::IParallelManager& pm) = 0;
};

} // namespace num
