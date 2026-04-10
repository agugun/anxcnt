#pragma once
#include "lib/modules.hpp"

namespace mod {
using namespace top;
namespace physics_heat {

class Heat3DState : public IState {
public:
    Vector temperatures;
    int nx, ny, nz;
    double dx, dy, dz;

    Heat3DState(int nx, int ny, int nz, double dx, double dy, double dz, double initial_temp = 0.0) 
        : nx(nx), ny(ny), nz(nz), dx(dx), dy(dy), dz(dz) {
        temperatures.resize(nx * ny * nz, initial_temp);
    }

    void update(const Vector& delta) override {
        for (size_t i = 0; i < temperatures.size(); ++i) {
            temperatures[i] += delta[i];
        }
    }

    std::unique_ptr<IState> clone() const override {
        auto copy = std::make_unique<Heat3DState>(nx, ny, nz, dx, dy, dz, 0.0);
        copy->temperatures = this->temperatures;
        return copy;
    }

    // Helper to get index
    int idx(int i, int j, int k) const {
        return (k * ny + j) * nx + i;
    }
};

} // namespace physics_heat
} // namespace mod
