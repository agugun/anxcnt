#pragma once
#include "lib/base.hpp"

namespace numerical_methods {
namespace physics_heat {

class Heat2DState : public IState {
public:
    Vector temperatures;
    int nx, ny;
    double dx, dy;

    Heat2DState(int nx, int ny, double dx, double dy, double initial_temp = 0.0) 
        : nx(nx), ny(ny), dx(dx), dy(dy) {
        temperatures.resize(nx * ny, initial_temp);
    }

    void update(const Vector& delta) override {
        for (size_t i = 0; i < temperatures.size(); ++i) {
            temperatures[i] += delta[i];
        }
    }

    std::unique_ptr<IState> clone() const override {
        auto copy = std::make_unique<Heat2DState>(nx, ny, dx, dy, 0.0);
        copy->temperatures = this->temperatures;
        return copy;
    }

    // Helper to get index
    int idx(int i, int j) const {
        return j * nx + i;
    }
};

} // namespace physics_heat
} // namespace numerical_methods
