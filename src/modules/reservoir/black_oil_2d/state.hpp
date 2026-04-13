#pragma once
#include "lib/spatial.hpp"
#include "lib/modules.hpp"
#include <vector>
#include <algorithm>

namespace mod {
using namespace top;
namespace reservoir {

/**
 * @brief ReservoirBlackOil2DState tracks 3 variables per cell.
 * Primary Variables: Oil Pressure [psi], Water Saturation [fraction], Gas Saturation [fraction].
 */
class ReservoirBlackOil2DState : public IState {
public:
    Spatial2D spatial;
    
    // Variables packed in a single vector [P, Sw, Sg, P, Sw, Sg, ...]
    Vector variables;

    ReservoirBlackOil2DState(Spatial2D spatial)
        : spatial(spatial),
          variables(3 * spatial.nx * spatial.ny, 0.0) {
        
        for (int i = 0; i < spatial.nx * spatial.ny; ++i) {
            variables[3 * i]     = 2000.0; // Initial Pressure
            variables[3 * i + 1] = 0.2;    // Initial Water (Connate)
            variables[3 * i + 2] = 0.1;    // Initial Gas
        }
    }

    // Accessors
    double p(int c) const { return variables[3 * c]; }
    double sw(int c) const { return variables[3 * c + 1]; }
    double sg(int c) const { return variables[3 * c + 2]; }
    double so(int c) const { return 1.0 - sw(c) - sg(c); }

    int idx(int i, int j) const { return spatial.idx(i, j); }

    void update(const Vector& delta) override {
        for (size_t i = 0; i < variables.size(); i += 3) {
            variables[i]     += delta[i];
            variables[i + 1] += delta[i + 1];
            variables[i + 2] += delta[i + 2];
            
            // Physical boundaries
            variables[i]     = std::max(14.7, variables[i]);
            variables[i + 1] = std::max(0.0, std::min(1.0, variables[i + 1]));
            variables[i + 2] = std::max(0.0, std::min(1.0 - variables[i + 1], variables[i + 2]));
        }
    }

    std::unique_ptr<IState> clone() const override {
        return std::make_unique<ReservoirBlackOil2DState>(*this);
    }
};

} // namespace reservoir
} // namespace mod
