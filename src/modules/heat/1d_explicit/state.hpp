#pragma once
#include "lib/spatial.hpp"
#include "../heat_state.hpp"

namespace mod {
using namespace top;
namespace heat {

class Heat1DExplicitState : public HeatState {
public:
    Spatial1D spatial;

    Heat1DExplicitState(Spatial1D spatial, double initial_temp) 
        : HeatState(spatial.nx, initial_temp), spatial(spatial) {}

    std::unique_ptr<IState> clone() const override {
        auto copy = std::make_unique<Heat1DExplicitState>(spatial, 0.0);
        copy->temperatures = this->temperatures;
        return copy;
    }
};

} // namespace heat
} // namespace mod