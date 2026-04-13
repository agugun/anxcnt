/**
 * @file math_utils.hpp
 * @brief Auxiliary Toolbox for Specialized Algorithms.
 * 
 * Objective:
 * This file provides specialized mathematical functions (FFT, Power Iteration) 
 * that are supplementary to the core simulation loop.
 * 
 * Architectural Rationale:
 * Segregating these utilities isolates non-core algorithmic complexities 
 * from the simulation's architectural pipeline. This keeps the primary 
 * library lightweight and prevents the main framework from becoming 
 * cluttered with algorithms that are only required for specialized 
 * analysis or diagnostics.
 * 
 * Strategic Importance:
 * By isolating these algorithms, we ensure that the core framework doesn't 
 * incur the maintenance cost of specialized numerical methods unless 
 * they are explicitly needed. It protects the simulation's "hot-path" from 
 * unnecessary dependency bloat.
 */
#pragma once
#include <vector>
#include <complex>
#include <cmath>
#include <functional>
#include <algorithm>
#include <stdexcept>

namespace num {

using Vector = std::vector<double>;
using VectorC = std::vector<std::complex<double>>;

/**
 * @brief Estimates the spectral radius (largest eigenvalue magnitude) using Power Iteration.
 * Useful for automated CFL stability analysis.
 * 
 * @param apply_A Function that performs matrix-vector multiplication (A * v).
 * @param n Size of the vector.
 * @param max_iter Maximum number of iterations (default: 100).
 * @param tol Convergence tolerance for the eigenvalue (default: 1e-4).
 * @return The estimated spectral radius.
 */
inline double estimate_spectral_radius(std::function<Vector(const Vector&)> apply_A, 
                                       size_t n, 
                                       int max_iter = 100, 
                                       double tol = 1e-4) {
    if (n == 0) return 0.0;

    // Start with a random vector
    Vector v(n);
    for (size_t i = 0; i < n; ++i) v[i] = (double)rand() / RAND_MAX;

    double lambda_old = 0.0;
    for (int i = 0; i < max_iter; ++i) {
        // Normalize
        double norm = 0.0;
        for (double val : v) norm += val * val;
        norm = std::sqrt(norm);
        if (norm < 1e-15) break; 
        for (double& val : v) val /= norm;

        // Apply operator
        Vector w = apply_A(v);

        // Rayleigh Quotient / Eigenvalue estimate: lambda = (v*w) / (v*v)
        // Since v is normalized, lambda = v*w
        double lambda = 0.0;
        for (size_t j = 0; j < n; ++j) lambda += v[j] * w[j];

        if (std::abs(lambda - lambda_old) < tol * std::abs(lambda)) {
            return std::abs(lambda);
        }

        v = w;
        lambda_old = lambda;
    }

    return std::abs(lambda_old);
}

/**
 * @brief Recursive Cooley-Tukey FFT implementation.
 * Requirements: Data size must be a power of two.
 */
inline void fft_recursive(VectorC& x, bool inverse) {
    size_t n = x.size();
    if (n <= 1) return;

    // Split into even and odd
    VectorC even(n / 2), odd(n / 2);
    for (size_t i = 0; i < n / 2; ++i) {
        even[i] = x[2 * i];
        odd[i] = x[2 * i + 1];
    }

    // Recursive calls
    fft_recursive(even, inverse);
    fft_recursive(odd, inverse);

    // Combine
    double angle = 2.0 * M_PI / n * (inverse ? 1 : -1);
    std::complex<double> w(1), wn(std::cos(angle), std::sin(angle));
    for (size_t i = 0; i < n / 2; ++i) {
        std::complex<double> t = w * odd[i];
        x[i] = even[i] + t;
        x[i + n / 2] = even[i] - t;
        w *= wn;
    }
}

/**
 * @brief Computes the Fast Fourier Transform of a complex vector.
 * If size is not a power of two, it will be padded with zeros internally.
 */
inline VectorC fft(const VectorC& input) {
    size_t n = input.size();
    if (n == 0) return {};

    // Pad to next power of two
    size_t m = 1;
    while (m < n) m <<= 1;
    
    VectorC x(m, 0.0);
    std::copy(input.begin(), input.end(), x.begin());

    fft_recursive(x, false);
    return x;
}

/**
 * @brief Computes the Inverse Fast Fourier Transform of a complex vector.
 */
inline VectorC ifft(const VectorC& input) {
    size_t n = input.size();
    if (n == 0) return {};

    VectorC x = input;
    fft_recursive(x, true);

    // Normalize
    for (auto& val : x) val /= (double)n;
    return x;
}

} // namespace num
