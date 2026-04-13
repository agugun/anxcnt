#pragma once
#include "lib/spatial.hpp"
#include "lib/modules.hpp"
#include <vector>
#include <memory>
#include <algorithm>

namespace mod {
using namespace top;
namespace reservoir {

class ReservoirOilGas2DState : public IState {
public:
    Spatial2D spatial;
    std::vector<double> pressures;   // Oil Phase Pressure [psi]
    std::vector<double> gas_saturations; // Gas Saturation [fraction]
    
    // Fixed Connate Water
    const double swc = 0.2;

    ReservoirOilGas2DState(Spatial2D spatial, double p_init, double s_init)
        : spatial(spatial), 
          pressures(spatial.nx * spatial.ny, p_init), 
          gas_saturations(spatial.nx * spatial.ny, s_init) {}

    int idx(int i, int j) const { return j * spatial.nx + i; }

    void update(const Vector& delta) override {
        for (size_t i = 0; i < pressures.size(); ++i) {
            pressures[i] += delta[2 * i];
            gas_saturations[i] += delta[2 * i + 1];
            
            // Physical boundaries
            pressures[i] = std::max(14.7, pressures[i]); // Cap at atmospheric
            gas_saturations[i] = std::max(0.0, std::min(1.0 - swc, gas_saturations[i]));
        }
    }

    std::unique_ptr<IState> clone() const override {
        return std::make_unique<ReservoirOilGas2DState>(*this);
    }
};

} // namespace reservoir
} // namespace mod
