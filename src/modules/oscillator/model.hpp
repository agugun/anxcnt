#pragma once
#include "lib/interfaces.hpp"
#include "state.hpp"

namespace mod::oscillator {

// Strong semantic typing
using Mass_kg = double;
using Damping_Ns_m = double;
using Stiffness_N_m = double;

/**
 * @brief Physics Model for the Harmonic Oscillator.
 * Only contains continuous physical properties and weight mappings.
 */
class OscillatorModel : public top::IModel {
public:
    Mass_kg m;
    Damping_Ns_m c;
    Stiffness_N_m k;

    OscillatorModel(Mass_kg mass, Damping_Ns_m damp, Stiffness_N_m stiff) 
        : m(mass), c(damp), k(stiff) {}

    double get_tolerance() const override { return 1e-6; }

    top::Vector get_accumulation_weights(const top::IGrid& grid, const top::IState& state) const override {
        // [Eq. 1 - Notebook] dx/dt = v  -> Weight = 1.0
        // [Eq. 2 - Notebook] m*dv/dt = ... -> Weight = m
        return {1.0, m};
    }
};

} // namespace mod::oscillator
