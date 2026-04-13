#pragma once
#include "lib/spatial.hpp"
#include "lib/modules.hpp"
#include <vector>

namespace mod {
using namespace top;
namespace reservoir {

class Reservoir1DState : public IState {
public:
    std::vector<double> pressures; // [psi]
    Spatial1D spatial;

    Reservoir1DState(Spatial1D spatial, double initial_p)
        : spatial(spatial), pressures(spatial.nx, initial_p) {}

    size_t size() const { return pressures.size(); }
    
    void update(const Vector& delta) override {
        if (delta.size() != pressures.size()) {
            throw std::runtime_error("State Update Error: Size mismatch between delta and pressures.");
        }
        for (size_t i = 0; i < pressures.size(); ++i) {
            pressures[i] += delta[i];
        }
    }

    std::unique_ptr<IState> clone() const override {
        auto copy = std::make_unique<Reservoir1DState>(spatial, 0.0);
        copy->pressures = this->pressures;
        return copy;
    }
};

} // namespace reservoir
} // namespace mod
