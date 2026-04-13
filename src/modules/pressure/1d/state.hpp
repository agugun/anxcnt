#pragma once
#include "lib/spatial.hpp"
#include "lib/modules.hpp"

namespace mod {
using namespace top;
namespace pressure {

class Pressure1DState : public IState {
public:
    Vector pressures;
    Spatial1D spatial;

    Pressure1DState(Spatial1D spatial, double initial_press) 
        : pressures(spatial.nx, initial_press), spatial(spatial) {}

    void update(const Vector& delta) override {
        for (size_t i = 0; i < pressures.size(); ++i) {
            pressures[i] += delta[i];
        }
    }

    std::unique_ptr<IState> clone() const override {
        auto copy = std::make_unique<Pressure1DState>(spatial, 0.0);
        copy->pressures = this->pressures;
        return copy;
    }
};

} // namespace pressure
} // namespace mod
