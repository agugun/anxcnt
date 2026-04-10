#pragma once
#include "lib/modules.hpp"

namespace mod {
using namespace top;
namespace physics_heat {

class Heat1DState : public IState {
public:
    Vector temperatures;
    double dx;

    Heat1DState(size_t size, double dx, double initial_temp) 
        : temperatures(size, initial_temp), dx(dx) {}

    void update(const Vector& delta) override {
        // Assume delta corresponds to interior points or full vector?
        // To follow base.hpp pattern exactly, delta should be the same size as temperatures.
        // However, in implicit solvers, often we only solve for interior.
        // I'll stick to full vector delta for consistency with IState::update(const Vector& delta).
        for (size_t i = 0; i < temperatures.size(); ++i) {
            temperatures[i] += delta[i];
        }
    }

    std::unique_ptr<IState> clone() const override {
        auto copy = std::make_unique<Heat1DState>(temperatures.size(), dx, 0.0);
        copy->temperatures = this->temperatures;
        return copy;
    }
};

} // namespace physics_heat
} // namespace mod
