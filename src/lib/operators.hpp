/**
 * @file operators.hpp
 * @brief Mathematical Toolbox for Discrete Calculus.
 * 
 * Objective:
 * This file centralizes the implementations of common numerical operators 
 * (Laplacian, Gradient) used across multiple physics modules.
 * 
 * Architectural Rationale:
 * Segregating discrete operators into a "Numerical Toolbox" ensures 
 * mathematical consistency. It prevents "Stencil Drift"—where different 
 * modules might implement slightly different version's of the same physics—
 * and allows for centralized performance optimization of numerical kernels.
 * 
 * Strategic Importance:
 * By providing a unified interface for spatial derivatives and external forces 
 * (Sources/Sinks), this file allows physics modules to focus purely on 
 * governing laws ($du/dt = ...$) while outsourcing the mechanics of 
 * spatial discretization.
 */
#pragma once
#include <vector>
#include <functional>
#include "lib/modules.hpp"

namespace top {

/**
 * @brief ISourceSink interface for external terms (sinks/sources).
 * Provides a standard way to apply external physics to residual and jacobian.
 */
class ISourceSink {
public:
    virtual ~ISourceSink() = default;

    /**
     * @brief Apply source term to the residual and jacobian.
     * @param residual The residual vector to modify (length 1*N or 2*N)
     * @param jacobian Pointer to jacobian matrix (can be nullptr)
     * @param state Current state of the field
     * @param dt Time step
     */
    virtual void apply(Vector& residual, Matrix* jacobian, const top::IState& state, double dt) = 0;
};

} // namespace top

namespace mop {
using namespace top;


using Vector = std::vector<double>;

/**
 * @brief 1D Central Difference Laplacian: (u[i+1] - 2u[i] + u[i-1]) / dx^2
 * 
 * Boundary: Simple Dirichlet assumption (zero at boundaries).
 * Complexity: O(N)
 */
inline Vector laplace_1d(const Vector& u, double dx) {
    size_t n = u.size();
    Vector res(n, 0.0);
    double inv_dx2 = 1.0 / (dx * dx);

    for (size_t i = 1; i < n - 1; ++i) {
        res[i] = (u[i+1] - 2.0 * u[i] + u[i-1]) * inv_dx2;
    }
    return res;
}

/**
 * @brief 2D 5-point Stencil Laplacian.
 * 
 * Discretization: d2u/dx2 + d2u/dy2
 * Fit With: Uniform Cartesian Grids.
 */
inline Vector laplace_2d(const Vector& u, int nx, int ny, double dx, double dy) {
    Vector res(u.size(), 0.0);
    double inv_dx2 = 1.0 / (dx * dx);
    double inv_dy2 = 1.0 / (dy * dy);

    auto idx = [nx](int i, int j) { return j * nx + i; };

    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            double d2u_dx2 = (u[idx(i + 1, j)] - 2.0 * u[idx(i, j)] + u[idx(i - 1, j)]) * inv_dx2;
            double d2u_dy2 = (u[idx(i, j + 1)] - 2.0 * u[idx(i, j)] + u[idx(i, j - 1)]) * inv_dy2;
            res[idx(i, j)] = d2u_dx2 + d2u_dy2;
        }
    }
    return res;
}

/**
 * @brief 1D Central Gradient: (u[i+1] - u[i-1]) / (2dx)
 */
inline Vector grad_1d(const Vector& u, double dx) {
    size_t n = u.size();
    Vector res(n, 0.0);
    double inv_2dx = 1.0 / (2.0 * dx);

    for (size_t i = 1; i < n - 1; ++i) {
        res[i] = (u[i+1] - u[i-1]) * inv_2dx;
    }
    return res;
}

/**
 * @brief 3D 7-point Stencil Laplacian.
 * 
 * Order of Accuracy: 2nd Order in Space.
 * Boundary: Zero boundary condition assumed for stencil neighbors.
 */
inline Vector laplace_3d(const Vector& u, int nx, int ny, int nz, double dx, double dy, double dz) {
    Vector res(u.size(), 0.0);
    double inv_dx2 = 1.0 / (dx * dx);
    double inv_dy2 = 1.0 / (dy * dy);
    double inv_dz2 = 1.0 / (dz * dz);

    auto idx = [nx, ny](int i, int j, int k) { 
        return (k * ny + j) * nx + i; 
    };

    for (int k = 1; k < nz - 1; ++k) {
        for (int j = 1; j < ny - 1; ++j) {
            for (int i = 1; i < nx - 1; ++i) {
                double d2u_dx2 = (u[idx(i + 1, j, k)] - 2.0 * u[idx(i, j, k)] + u[idx(i - 1, j, k)]) * inv_dx2;
                double d2u_dy2 = (u[idx(i, j + 1, k)] - 2.0 * u[idx(i, j, k)] + u[idx(i, j - 1, k)]) * inv_dy2;
                double d2u_dz2 = (u[idx(i, j, k + 1)] - 2.0 * u[idx(i, j, k)] + u[idx(i, j, k - 1)]) * inv_dz2;
                res[idx(i, j, k)] = d2u_dx2 + d2u_dy2 + d2u_dz2;
            }
        }
    }
    return res;
}

} // namespace mop
