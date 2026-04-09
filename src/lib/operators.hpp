#pragma once
#include <vector>
#include <functional>

namespace numerical_methods {
namespace operators {

using Vector = std::vector<double>;

// 1D Central Difference Laplacian: (u[i+1] - 2u[i] + u[i-1]) / dx^2
inline Vector laplace_1d(const Vector& u, double dx) {
    size_t n = u.size();
    Vector res(n, 0.0);
    double inv_dx2 = 1.0 / (dx * dx);

    for (size_t i = 1; i < n - 1; ++i) {
        res[i] = (u[i+1] - 2.0 * u[i] + u[i-1]) * inv_dx2;
    }
    return res;
}

// 2D 5-point Stencil Laplacian
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

// 1D Gradient (Central Difference): (u[i+1] - u[i-1]) / (2dx)
inline Vector grad_1d(const Vector& u, double dx) {
    size_t n = u.size();
    Vector res(n, 0.0);
    double inv_2dx = 1.0 / (2.0 * dx);

    for (size_t i = 1; i < n - 1; ++i) {
        res[i] = (u[i+1] - u[i-1]) * inv_2dx;
    }
    return res;
}

} // namespace operators
} // namespace numerical_methods
