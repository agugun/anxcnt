/**
 * @file iterative_solvers.hpp
 * @brief Scalable Iterative Linear Solvers for Sparse Systems.
 */
#pragma once
#include "sparse.hpp"
#include <vector>
#include <cmath>
#include <functional>
#include <iostream>

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
        
        // Memory Buffers (Allocate once)
        std::vector<double> r(n), r_hat(n), v(n), p(n), s(n), t(n);
        
        // initial r = b - A*x
        std::vector<double> Ax0 = A.multiply(x);
        #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            r[i] = b[i] - Ax0[i];
            r_hat[i] = r[i];
            v[i] = 0.0;
            p[i] = 0.0;
        }
        
        double rho = 1.0, alpha = 1.0, omega = 1.0;
        
        for (int i = 0; i < max_iter; ++i) {
            double rho_new = 0.0;
            #pragma omp parallel for reduction(+:rho_new)
            for(int j=0; j<n; ++j) rho_new += r_hat[j] * r[j];
            
            if (rho_new != rho_new) { // Check NaN
                if (verbose) std::cerr << "  [BiCGSTAB] Error: rho_new is NaN" << std::endl;
                break;
            }

            if (std::abs(rho_new) < 1e-60) {
                if (verbose) std::cerr << "  [BiCGSTAB] Warning: rho_new too small: " << rho_new << std::endl;
                break;
            }
            
            if (i == 0) {
                #pragma omp parallel for
                for(int j=0; j<n; ++j) p[j] = r[j];
            } else {
                double beta = (rho_new / rho) * (alpha / omega);
                if (std::isnan(beta) || std::isinf(beta)) {
                   if (verbose) std::cerr << "  [BiCGSTAB] Error: beta is NaN/Inf at iter " << i << ". rho_new=" << rho_new << " rho=" << rho << " alpha=" << alpha << " omega=" << omega << std::endl;
                   break;
                }
                #pragma omp parallel for
                for(int j=0; j<n; ++j) {
                    p[j] = r[j] + beta * (p[j] - omega * v[j]);
                }
            }
            
            v = A.multiply(p);
            double v_dot_rhat = 0.0;
            #pragma omp parallel for reduction(+:v_dot_rhat)
            for(int j=0; j<n; ++j) v_dot_rhat += r_hat[j] * v[j];
            
            if (std::abs(v_dot_rhat) < 1e-60) {
                if (verbose) std::cerr << "  [BiCGSTAB] Warning: v_dot_rhat too small: " << v_dot_rhat << std::endl;
                break;
            }

            alpha = rho_new / v_dot_rhat;
            
            double norm_s = 0.0;
            #pragma omp parallel for reduction(+:norm_s)
            for(int j=0; j<n; ++j) {
                s[j] = r[j] - alpha * v[j];
                norm_s += s[j] * s[j];
            }
            norm_s = std::sqrt(norm_s);
            
            if (verbose && std::isnan(norm_s)) {
                std::cerr << "  [BiCGSTAB] Error: norm_s is NaN. rho_new=" << rho_new << " v_dot_rhat=" << v_dot_rhat << " alpha=" << alpha << std::endl;
                break;
            }
            
            if (norm_s < tol) {
                #pragma omp parallel for
                for(int j=0; j<n; ++j) x[j] += alpha * p[j];
                if (verbose) std::cout << "  [BiCGSTAB] Conv after " << i << " iterations. Res=" << norm_s << std::endl;
                return x;
            }
            
            t = A.multiply(s);
            double t_dot_t = 0.0, t_dot_s = 0.0;
            #pragma omp parallel for reduction(+:t_dot_t, t_dot_s)
            for(int j=0; j<n; ++j) {
                t_dot_t += t[j] * t[j];
                t_dot_s += t[j] * s[j];
            }
            
            omega = (t_dot_t < 1e-60) ? 0.0 : t_dot_s / t_dot_t;
            
            double norm_r = 0.0;
            #pragma omp parallel for reduction(+:norm_r)
            for(int j=0; j<n; ++j) {
                x[j] += alpha * p[j] + omega * s[j];
                r[j] = s[j] - omega * t[j];
                norm_r += r[j] * r[j];
            }
            norm_r = std::sqrt(norm_r);
            
            if (norm_r < tol) {
                if (verbose) std::cout << "  [BiCGSTAB] Conv after " << i << " iterations. Res=" << norm_r << std::endl;
                return x;
            }
            rho = rho_new;
            if (verbose && i % 50 == 0) std::cout << "    L-Iter " << i << ": Res=" << norm_r << std::endl;
        }
        
        return x;
    }
};

} // namespace num
