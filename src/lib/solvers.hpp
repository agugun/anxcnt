#pragma once
#include "lib/base.hpp"
#include <stdexcept>
#include <functional>
#include <cmath>

namespace numerical_methods {

/**
 * @brief LinearTridiagonalSolver implements the ISolver interface for 1D problems.
 * It strictly uses the base.hpp framework's Matrix and Vector types.
 */
class LinearTridiagonalSolver : public ISolver {
public:
    Vector solve(const IModel& model, const IState& state, double dt) override {
        // This method is a placeholder for the generic ISolver interface.
        // For shared implicit logic, see solve_system below.
        return {};
    }

    /**
     * @brief Specialized Thomas Algorithm for tridiagonal systems.
     * Extracts J diagonals and solves J * delta = -res.
     */
    Vector solve_system(const Matrix& J, const Vector& res) {
        size_t n = res.size();
        if (n == 0) return {};

        Vector a(n, 0.0), b(n, 0.0), c(n, 0.0), d(n, 0.0);
        for (size_t i = 0; i < n; ++i) {
            b[i] = J[i][i];
            if (i > 0) a[i] = J[i][i-1];
            if (i < n - 1) c[i] = J[i][i+1];
            d[i] = -res[i];
        }

        Vector c_prime(n);
        Vector d_prime(n);
        Vector x(n);

        // Forward elimination
        c_prime[0] = c[0] / b[0];
        d_prime[0] = d[0] / b[0];

        for (size_t i = 1; i < n; ++i) {
            double den = b[i] - a[i] * c_prime[i-1];
            if (den == 0) throw std::runtime_error("TDMA: Divide by zero");
            double m = 1.0 / den;
            c_prime[i] = c[i] * m;
            d_prime[i] = (d[i] - a[i] * d_prime[i-1]) * m;
        }

        // Backward substitution
        x[n-1] = d_prime[n-1];
        for (int i = (int)n - 2; i >= 0; --i) {
            x[i] = d_prime[i] - c_prime[i] * x[i + 1];
        }

        return x;
    }
};

/**
 * @brief ConjugateGradientSolver for symmetric positive definite systems.
 * Useful for 2D/3D implicit problems where dense matrices are too large.
 */
class ConjugateGradientSolver : public ISolver {
public:
    Vector solve(const IModel& model, const IState& state, double dt) override {
        return {}; // Placeholder
    }

    /**
     * @brief Solves A*x = b using the Conjugate Gradient method.
     * 'apply_A' is a callback that calculates the matrix-vector product.
     */
    Vector solve_iterative(std::function<Vector(const Vector&)> apply_A, 
                           const Vector& b, 
                           const Vector& x0, 
                           double tol = 1e-6, 
                           int max_iter = 1000) {
        Vector x = x0;
        Vector r = b;
        Vector Ax = apply_A(x);
        for(size_t i=0; i<r.size(); ++i) r[i] -= Ax[i];
        
        Vector p = r;
        double rsold = dot(r, r);

        for (int i = 0; i < max_iter; ++i) {
            Vector Ap = apply_A(p);
            double alpha = rsold / dot(p, Ap);
            
            for(size_t j=0; j<x.size(); ++j) {
                x[j] += alpha * p[j];
                r[j] -= alpha * Ap[j];
            }
            
            double rsnew = dot(r, r);
            if (std::sqrt(rsnew) < tol) {
                // std::cout << "CG converged in " << i << " iterations.\n";
                break;
            }
            
            for(size_t j=0; j<p.size(); ++j) {
                p[j] = r[j] + (rsnew / rsold) * p[j];
            }
            rsold = rsnew;
        }
        return x;
    }

private:
    double dot(const Vector& v1, const Vector& v2) {
        double res = 0;
        for(size_t i=0; i<v1.size(); ++i) res += v1[i] * v2[i];
        return res;
    }
};

} // namespace numerical_methods
