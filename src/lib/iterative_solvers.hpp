/**
 * @file iterative_solvers.hpp
 * @brief Scalable Iterative Linear Solvers for Sparse Systems.
 */
#pragma once
#include "sparse.hpp"
#include <vector>
#include <cmath>
#include <functional>

namespace num {

/**
 * @brief BiCGSTAB (Biconjugate Gradient Stabilized) Solver.
 * 
 * Fits With:
 * - Large, non-symmetric sparse systems.
 * - Non-linear reservoir models where assembling full LU is prohibitive.
 */
class BiCGSTABSolver {
public:
    static std::vector<double> solve(const SparseMatrix& A, 
                                     const std::vector<double>& b, 
                                     double tol = 1e-8, 
                                     int max_iter = 1000,
                                     bool verbose = false) {
        int n = (int)b.size();
        std::vector<double> x(n, 0.0);
        std::vector<double> r = subtract(b, A.multiply(x));
        std::vector<double> r_hat = r;
        
        double rho = 1.0, alpha = 1.0, omega = 1.0;
        std::vector<double> v(n, 0.0), p(n, 0.0);
        
        for (int i = 0; i < max_iter; ++i) {
            double rho_new = dot(r_hat, r);
            if (std::abs(rho_new) < 1e-20) break;
            
            if (i == 0) {
                p = r;
            } else {
                double beta = (rho_new / rho) * (alpha / omega);
                p = add(r, scale(subtract(p, scale(v, omega)), beta));
            }
            
            v = A.multiply(p);
            double v_dot_rhat = dot(r_hat, v);
            if (std::abs(v_dot_rhat) < 1e-20) break;

            alpha = rho_new / v_dot_rhat;
            std::vector<double> s = subtract(r, scale(v, alpha));
            
            double norm_s = norm(s);
            if (norm_s < tol) {
                x = add(x, scale(p, alpha));
                if (verbose) std::cout << "  [BiCGSTAB] Conv after " << i << " iterations. Res=" << norm_s << std::endl;
                return x;
            }
            
            std::vector<double> t = A.multiply(s);
            double t_dot_t = dot(t, t);
            omega = (t_dot_t < 1e-20) ? 0.0 : dot(t, s) / t_dot_t;
            
            x = add(x, add(scale(p, alpha), scale(s, omega)));
            r = subtract(s, scale(t, omega));
            
            double norm_r = norm(r);
            if (norm_r < tol) {
                if (verbose) std::cout << "  [BiCGSTAB] Conv after " << i << " iterations. Res=" << norm_r << std::endl;
                return x;
            }
            rho = rho_new;
            if (verbose && i % 50 == 0) std::cout << "    L-Iter " << i << ": Res=" << norm_r << std::endl;
        }
        
        return x;
    }

private:
    static double dot(const std::vector<double>& a, const std::vector<double>& b) {
        double res = 0;
        for (size_t i = 0; i < a.size(); ++i) res += a[i] * b[i];
        return res;
    }
    
    static double norm(const std::vector<double>& v) {
        return std::sqrt(dot(v, v));
    }
    
    static std::vector<double> add(std::vector<double> a, const std::vector<double>& b) {
        for (size_t i = 0; i < a.size(); ++i) a[i] += b[i];
        return a;
    }
    
    static std::vector<double> subtract(std::vector<double> a, const std::vector<double>& b) {
        for (size_t i = 0; i < a.size(); ++i) a[i] -= b[i];
        return a;
    }
    
    static std::vector<double> scale(std::vector<double> a, double s) {
        for (size_t i = 0; i < a.size(); ++i) a[i] *= s;
        return a;
    }
};

} // namespace num
