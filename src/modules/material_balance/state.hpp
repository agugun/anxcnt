#pragma once
#include "lib/base.hpp"

namespace numerical_methods {
namespace physics_mba {

class MBState : public IState {
public:
    double pressure;

    MBState(double initial_pressure) : pressure(initial_pressure) {}

    void update(const Vector& delta) override {
        // Delta is a vector of size 1
        pressure += delta[0];
    }

    std::unique_ptr<IState> clone() const override {
        return std::make_unique<MBState>(pressure);
    }
};

} // namespace physics_mba
} // namespace numerical_methods
