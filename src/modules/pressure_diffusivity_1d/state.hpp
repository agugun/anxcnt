#pragma once
#include "lib/base.hpp"

namespace numerical_methods {
namespace physics_pressure {

class Pressure1DState : public IState {
public:
    Vector pressures;
    double dx;

    Pressure1DState(size_t size, double dx, double initial_press) 
        : pressures(size, initial_press), dx(dx) {}

    void update(const Vector& delta) override {
        for (size_t i = 0; i < pressures.size(); ++i) {
            pressures[i] += delta[i];
        }
    }

    std::unique_ptr<IState> clone() const override {
        auto copy = std::make_unique<Pressure1DState>(pressures.size(), dx, 0.0);
        copy->pressures = this->pressures;
        return copy;
    }
};

} // namespace physics_pressure
} // namespace numerical_methods
