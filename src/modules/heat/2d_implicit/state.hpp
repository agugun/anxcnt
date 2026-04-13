#pragma once
#include "lib/spatial.hpp"
#include "../heat_state.hpp"

namespace mod {
using namespace top;
namespace heat {

class Heat2DImplicitState : public HeatState {
public:
    Spatial2D spatial;

    Heat2DImplicitState(Spatial2D spatial, double initial_temp = 0.0) 
        : HeatState(spatial.total_size(), initial_temp), spatial(spatial) {}

    std::unique_ptr<IState> clone() const override {
        auto copy = std::make_unique<Heat2DImplicitState>(spatial, 0.0);
        copy->temperatures = this->temperatures;
        return copy;
    }
};

} // namespace heat
} // namespace mod
