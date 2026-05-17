#pragma once
#include "lib/interfaces.hpp"

namespace mod::reservoir {

/**
 * @brief Single-cell grid for 0D tank simulations.
 */
class MBAGrid : public IGrid {
public:
    size_t get_total_cells() const override { return 1; }
};

/**
 * @brief State representing a 0D Tank (Pressure only).
 */
class MBAState : public IState {
public:
    double pressure;

    MBAState(double initial_pressure) : pressure(initial_pressure) {}

    void update(const std::vector<double>& delta) override {
        pressure += delta[0];
    }

    std::vector<double> to_vector() const override {
        return { pressure };
    }

    std::unique_ptr<IState> clone() const override {
        return std::make_unique<MBAState>(pressure);
    }
};

} // namespace mod::reservoir
