#pragma once
#include "lib/spatial.hpp"
#include "../heat_state.hpp"

namespace mod {
using namespace top;
namespace heat {

class Heat3DImplicitState : public HeatState {
public:
    std::shared_ptr<Spatial3D> spatial;

    Heat3DImplicitState(std::shared_ptr<Spatial3D> s, double initial_temp = 0.0) 
        : HeatState(s->total_size(), initial_temp), spatial(s) {}

    std::unique_ptr<IState> clone() const override {
        auto copy = std::make_unique<Heat3DImplicitState>(spatial, 0.0);
        copy->temperatures = this->temperatures;
        return copy;
    }
};

} // namespace heat
} // namespace mod
