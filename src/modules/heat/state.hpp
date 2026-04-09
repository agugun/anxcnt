#pragma once
#include "lib/base.hpp"

class HeatState : public IState {
public:
    Vector temperatures;

    HeatState(size_t size, double initial_temp) : temperatures(size, initial_temp) {}

    void update(const Vector& delta) override {
        for (size_t i = 0; i < temperatures.size(); ++i) {
            temperatures[i] += delta[i];
        }
    }

    std::unique_ptr<IState> clone() const override {
        auto copy = std::make_unique<HeatState>(temperatures.size(), 0.0);
        copy->temperatures = this->temperatures;
        return copy;
    }
};