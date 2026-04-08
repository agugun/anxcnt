#ifndef TRIDIAGONAL_SOLVER_HPP
#define TRIDIAGONAL_SOLVER_HPP

#include <vector>
#include <stdexcept>

namespace numerical_methods {
namespace solver {

/**
 * @brief Thomas Algorithm (TDMA) for tridiagonal matrix inversion.
 * Solves: a[i]*x[i-1] + b[i]*x[i] + c[i]*x[i+1] = d[i]
 */
class TridiagonalSolver {
public:
    static std::vector<double> solve(const std::vector<double>& a, 
                                    const std::vector<double>& b, 
                                    const std::vector<double>& c, 
                                    const std::vector<double>& d) {
        int n = static_cast<int>(d.size());
        if (n == 0) return {};
        if (static_cast<int>(b.size()) != n || static_cast<int>(a.size()) != n || static_cast<int>(c.size()) != n) {
            throw std::invalid_argument("Vector sizes must match.");
        }

        std::vector<double> c_prime(n);
        std::vector<double> d_prime(n);
        std::vector<double> x(n);

        // Forward elimination
        c_prime[0] = c[0] / b[0];
        d_prime[0] = d[0] / b[0];

        for (int i = 1; i < n; ++i) {
            double m = 1.0 / (b[i] - a[i] * c_prime[i - 1]);
            c_prime[i] = c[i] * m;
            d_prime[i] = (d[i] - a[i] * d_prime[i - 1]) * m;
        }

        // Backward substitution
        x[n - 1] = d_prime[n - 1];
        for (int i = n - 2; i >= 0; --i) {
            x[i] = d_prime[i] - c_prime[i] * x[i + 1];
        }

        return x;
    }
};

} // namespace solver
} // namespace numerical_methods

#endif // TRIDIAGONAL_SOLVER_HPP
