#pragma once
#include <cmath>
#include <algorithm>

namespace mod {
namespace reservoir {

/**
 * @brief BlackOilPVT encapsulates industry-standard correlations (Standing's)
 * for Oil, Gas, and Water properties.
 */
class BlackOilPVT {
public:
    double gamma_g; // Gas specific gravity (air=1.0)
    double gamma_o; // Oil specific gravity (water=1.0)
    double temp_f;  // Reservoir Temperature [F]
    double pb;      // Bubble Point [psi]
    
    double api;     // API Gravity
    
    BlackOilPVT(double g_g = 0.7, double g_o = 0.8, double temp = 150.0, double p_b = 1500.0)
        : gamma_g(g_g), gamma_o(g_o), temp_f(temp), pb(p_b) {
        api = (141.5 / gamma_o) - 131.5;
    }

    // --- Solution Gas Oil Ratio (Rs) [SCF/STB] ---
    double get_rs(double p) const {
        double p_eff = std::min(pb, std::max(14.7, p));
        double a = 0.0125 * api - 0.00091 * temp_f;
        return gamma_g * std::pow((p_eff / 18.2) * std::pow(10, a), 1.2048);
    }

    double get_drs_dp(double p) const {
        if (p > pb) return 0.0;
        double p_capped = std::max(14.7, p);
        double a = 0.0125 * api - 0.00091 * temp_f;
        double coeff = gamma_g * std::pow(std::pow(10, a) / 18.2, 1.2048);
        return coeff * 1.2048 * std::pow(p_capped, 0.2048);
    }

    // --- Oil Formation Volume Factor (Bo) [rb/STB] ---
    double get_bo(double p, double rs) const {
        double f = rs * std::sqrt(gamma_g / gamma_o) + 1.25 * temp_f;
        double bo_saturated = 0.9759 + 0.00012 * std::pow(f, 1.2);
        if (p <= pb) return bo_saturated;
        return bo_saturated * std::exp(-1e-5 * (p - pb));
    }

    double get_dbo_dp(double p) const {
        double rs = get_rs(p);
        double f = rs * std::sqrt(gamma_g / gamma_o) + 1.25 * temp_f;
        double dbo_drs = 0.00012 * 1.2 * std::pow(f, 0.2) * std::sqrt(gamma_g / gamma_o);
        double sat_deriv = dbo_drs * get_drs_dp(p);
        if (p <= pb) return sat_deriv;
        return -1e-5 * get_bo(p, rs);
    }

    // --- Gas Formation Volume Factor (Bg) [rb/SCF] ---
    // Bg = 0.00503 * Z * T / P  (Assuming Z=1 for simplicity, or Z correlation)
    double get_bg(double p) const {
        double t_rankine = temp_f + 459.67;
        return 0.02827 * t_rankine / std::max(14.7, p); // Standard conversion factor
    }

    double get_dbg_dp(double p) const {
        double t_rankine = temp_f + 459.67;
        return -0.02827 * t_rankine / (std::max(14.7, p) * std::max(14.7, p));
    }

    // --- Water Formation Volume Factor (Bw) [rb/STB] ---
    double get_bw(double p) const {
        double cw = 3.0e-6; // Water compressibility
        return 1.0 / (1.0 + cw * (p - 14.7));
    }

    // --- Viscosities [cP] ---
    double get_mu_o(double p, double rs) const {
        // Simple dead oil viscosity + saturation adjustment
        double mu_dead = std::exp(std::exp(3.0 - 0.02 * api));
        return mu_dead / (1.0 + 0.001 * rs);
    }

    double get_mu_g(double p) const {
        return 0.013 + 2e-6 * p; // Simplified linear proxy
    }

    double get_mu_w(double p) const {
        return 0.5 + 1e-7 * p; // Simplified water viscosity
    }
};

} // namespace reservoir
} // namespace mod
