#pragma once
#include "lib/spatial.hpp"
#include "lib/modules.hpp"
#include <vector>
#include <memory>

namespace mod {
using namespace top;
namespace reservoir {

class ReservoirDualPhase2DState : public IState {
public:
    Spatial2D spatial;
    std::vector<double> pressures;         // Oil pressures [psi]
    std::vector<double> water_saturations; // Water saturations [fraction]

    ReservoirDualPhase2DState(Spatial2D spatial, double p_init, double s_init)
        : spatial(spatial),
          pressures(spatial.nx * spatial.ny, p_init),
          water_saturations(spatial.nx * spatial.ny, s_init) {}

    size_t size() const { return pressures.size(); }
    
    int idx(int i, int j) const { return j * spatial.nx + i; }

    // Flexible update for IMPES (pressure only) or Fully Implicit (coupled)
    void update(const Vector& delta) override {
        if (delta.size() == pressures.size()) {
            for (size_t i = 0; i < pressures.size(); ++i) {
                pressures[i] += delta[i];
            }
        } else if (delta.size() == 2 * pressures.size()) {
            for (size_t i = 0; i < pressures.size(); ++i) {
                pressures[i] += delta[2 * i];
                water_saturations[i] += delta[2 * i + 1];
                water_saturations[i] = std::max(0.0, std::min(1.0, water_saturations[i]));
            }
        }
    }

    std::unique_ptr<IState> clone() const override {
        auto copy = std::make_unique<ReservoirDualPhase2DState>(spatial, 0.0, 0.0);
        copy->pressures = this->pressures;
        copy->water_saturations = this->water_saturations;
        return copy;
    }
};

} // namespace reservoir
} // namespace mod
