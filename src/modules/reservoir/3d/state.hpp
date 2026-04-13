#pragma once
#include "lib/spatial.hpp"
#include "lib/modules.hpp"
#include <vector>
#include <memory>

namespace mod {
using namespace top;
namespace reservoir {

class Reservoir3DState : public IState {
public:
    std::vector<double> pressures; // [psi]
    Spatial3D spatial;

    Reservoir3DState(Spatial3D spatial, double initial_p)
        : spatial(spatial), pressures(spatial.nx * spatial.ny * spatial.nz, initial_p) {}

    size_t size() const { return pressures.size(); }
    
    int idx(int i, int j, int k) const { 
        return spatial.idx(i, j, k);
    }

    void update(const Vector& delta) override {
        for (size_t i = 0; i < pressures.size(); ++i) {
            pressures[i] += delta[i];
        }
    }

    std::unique_ptr<IState> clone() const override {
        auto copy = std::make_unique<Reservoir3DState>(spatial, 0.0);
        copy->pressures = this->pressures;
        return copy;
    }
};

} // namespace reservoir
} // namespace mod
