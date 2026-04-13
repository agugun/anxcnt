#pragma once
#include "lib/modules.hpp"

namespace mod {
using namespace top;
namespace heat {

/**
 * @brief Base class for all Heat Simulation states.
 * 
 * Centralizes the temperature vector and standard additive updates.
 */
class HeatState : public IState {
public:
    Vector temperatures;

    HeatState(size_t size, double initial_temp) : temperatures(size, initial_temp) {}
    
    virtual ~HeatState() = default;

    /**
     * @brief Update temperature values using a delta vector.
     */
    void update(const Vector& delta) override {
        if (delta.size() != temperatures.size()) {
            throw std::runtime_error("State Update Error: Size mismatch between delta and temperatures.");
        }
        for (size_t i = 0; i < temperatures.size(); ++i) {
            temperatures[i] += delta[i];
        }
    }

    // clone() remains pure virtual to be implemented by spatial subclasses
};

} // namespace heat
} // namespace mod
