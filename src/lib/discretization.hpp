#pragma once
#include <vector>
#include <memory>
#include <algorithm>
#include "interfaces.hpp"

namespace num {

/**
 * @namespace discretization
 * @brief Centralized Numerical Discretization Toolbox.
 */
namespace discretization {

/**
 * @brief 1D Conductance and Geometric properties.
 */
struct Conductance1D {
    std::vector<double> T; // Transmissibility
    Conductance1D(size_t n) : T(n - 1, 0.0) {}
};

/**
 * @brief 2D Conductance and Geometric properties.
 */
struct Conductance2D {
    std::vector<double> Tx, Ty;
    Conductance2D(size_t nx, size_t ny) : Tx((nx - 1) * ny, 0.0), Ty(nx * (ny - 1), 0.0) {}
};

/**
 * @brief 3D Conductance and Geometric properties.
 */
struct Conductance3D {
    std::vector<double> Tx, Ty, Tz;
    Conductance3D(size_t nx, size_t ny, size_t nz) 
        : Tx((nx - 1) * ny * nz, 0.0), Ty(nx * (ny - 1) * nz, 0.0), Tz(nx * ny * (nz - 1), 0.0) {}
};

// --- Heat Support (S.I. Units) ---

inline std::shared_ptr<Conductance1D> heat_cond_1d(int nx, double dx, double k, double h) {
    auto cond = std::make_shared<Conductance1D>(nx);
    double t_val = k * h / dx;
    std::fill(cond->T.begin(), cond->T.end(), t_val);
    return cond;
}

inline std::shared_ptr<Conductance2D> heat_cond_2d(int nx, int ny, double dx, double dy, double k, double thickness) {
    auto cond = std::make_shared<Conductance2D>(nx, ny);
    double tx_val = k * thickness * dy / dx;
    double ty_val = k * thickness * dx / dy;
    std::fill(cond->Tx.begin(), cond->Tx.end(), tx_val);
    std::fill(cond->Ty.begin(), cond->Ty.end(), ty_val);
    return cond;
}

inline std::shared_ptr<Conductance3D> heat_cond_3d(int nx, int ny, int nz, double dx, double dy, double dz, double k) {
    auto cond = std::make_shared<Conductance3D>(nx, ny, nz);
    double tx_val = k * dy * dz / dx;
    double ty_val = k * dx * dz / dy;
    double tz_val = k * dx * dy / dz;
    std::fill(cond->Tx.begin(), cond->Tx.end(), tx_val);
    std::fill(cond->Ty.begin(), cond->Ty.end(), ty_val);
    std::fill(cond->Tz.begin(), cond->Tz.end(), tz_val);
    return cond;
}

inline std::vector<double> heat_storage(int n, double vol, double rho, double cp) {
    return std::vector<double>(n, vol * rho * cp);
}

// --- Reservoir Support (Field Units) ---

inline std::shared_ptr<Conductance1D> reservoir_cond_1d(int nx, double dx, double k, double mu, double B, double h) {
    auto cond = std::make_shared<Conductance1D>(nx);
    double mult = 0.001127 * 5.615 * k / (mu * B);
    double t_val = mult * h / dx;
    std::fill(cond->T.begin(), cond->T.end(), t_val);
    return cond;
}

inline std::shared_ptr<Conductance2D> reservoir_cond_2d(int nx, int ny, double dx, double dy, double k, double mu, double B, double h) {
    auto cond = std::make_shared<Conductance2D>(nx, ny);
    double mult = 0.001127 * 5.615 * k / (mu * B);
    double tx_val = mult * dy * h / dx;
    double ty_val = mult * dx * h / dy;
    std::fill(cond->Tx.begin(), cond->Tx.end(), tx_val);
    std::fill(cond->Ty.begin(), cond->Ty.end(), ty_val);
    return cond;
}

inline std::shared_ptr<Conductance3D> reservoir_cond_3d(int nx, int ny, int nz, double dx, double dy, double dz, double k, double mu, double B) {
    auto cond = std::make_shared<Conductance3D>(nx, ny, nz);
    double mult = 0.001127 * 5.615 * k / (mu * B);
    double tx_val = mult * dy * dz / dx;
    double ty_val = mult * dx * dz / dy;
    double tz_val = mult * dx * dy / dz;
    std::fill(cond->Tx.begin(), cond->Tx.end(), tx_val);
    std::fill(cond->Ty.begin(), cond->Ty.end(), ty_val);
    std::fill(cond->Tz.begin(), cond->Tz.end(), tz_val);
    return cond;
}

inline std::vector<double> reservoir_storage(int n, double vol, double phi, double ct, double B) {
    return std::vector<double>(n, vol * phi * ct / B);
}

// --- Wave Support ---

inline std::shared_ptr<Conductance1D> wave_cond_1d(int nx, double dx, double c) {
    auto cond = std::make_shared<Conductance1D>(nx);
    double t_val = (c * c) / (dx * dx);
    std::fill(cond->T.begin(), cond->T.end(), t_val);
    return cond;
}

inline std::shared_ptr<Conductance2D> wave_cond_2d(int nx, int ny, double dx, double dy, double c) {
    auto cond = std::make_shared<Conductance2D>(nx, ny);
    double c2 = c * c;
    double tx_val = c2 / (dx * dx);
    double ty_val = c2 / (dy * dy);
    std::fill(cond->Tx.begin(), cond->Tx.end(), tx_val);
    std::fill(cond->Ty.begin(), cond->Ty.end(), ty_val);
    return cond;
}

inline std::vector<double> wave_storage(int n) {
    return std::vector<double>(n, 1.0);
}

// --- Pressure Support ---

inline std::shared_ptr<Conductance1D> pressure_cond_1d(int nx, double dx, double k, double mu, double area) {
    auto cond = std::make_shared<Conductance1D>(nx);
    double t_val = 0.001127 * k * area / (mu * dx);
    std::fill(cond->T.begin(), cond->T.end(), t_val);
    return cond;
}

inline std::shared_ptr<Conductance2D> pressure_cond_2d(int nx, int ny, double dx, double dy, double k, double mu, double area) {
    auto cond = std::make_shared<Conductance2D>(nx, ny);
    double t_val = 0.001127 * k * area / mu;
    double tx = t_val * dy / dx;
    double ty = t_val * dx / dy;
    std::fill(cond->Tx.begin(), cond->Tx.end(), tx);
    std::fill(cond->Ty.begin(), cond->Ty.end(), ty);
    return cond;
}

inline std::shared_ptr<Conductance3D> pressure_cond_3d(int nx, int ny, int nz, double dx, double dy, double dz, double k, double mu) {
    auto cond = std::make_shared<Conductance3D>(nx, ny, nz);
    double t_val = 0.001127 * k / mu;
    double tx = t_val * dy * dz / dx;
    double ty = t_val * dx * dz / dy;
    double tz = t_val * dx * dy / dz;
    std::fill(cond->Tx.begin(), cond->Tx.end(), tx);
    std::fill(cond->Ty.begin(), cond->Ty.end(), ty);
    std::fill(cond->Tz.begin(), cond->Tz.end(), tz);
    return cond;
}

inline std::vector<double> pressure_storage(int n, double vol, double phi, double ct) {
    return std::vector<double>(n, vol * phi * ct / 0.0002637);
}

} // namespace discretization
} // namespace num
