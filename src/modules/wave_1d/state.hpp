#pragma once
#include "lib/base.hpp"

namespace numerical_methods {
namespace physics_wave {

class Wave1DState : public IState {
public:
    Vector u; // Displacement
    Vector v; // Velocity
    double dx;

    Wave1DState(size_t size, double dx, double initial_displacement = 0.0) 
        : u(size, initial_displacement), v(size, 0.0), dx(dx) {}

    void update(const Vector& delta) override {
        // Assume delta is concatenated: [delta_u, delta_v]
        size_t n = u.size();
        for (size_t i = 0; i < n; ++i) {
            u[i] += delta[i];
            v[i] += delta[i + n];
        }
    }

    std::unique_ptr<IState> clone() const override {
        auto copy = std::make_unique<Wave1DState>(u.size(), dx, 0.0);
        copy->u = this->u;
        copy->v = this->v;
        return copy;
    }

    // Helper to get concatenated state for RHS evaluation
    Vector get_combined() const {
        size_t n = u.size();
        Vector combined(2 * n);
        for (size_t i = 0; i < n; ++i) {
            combined[i] = u[i];
            combined[i + n] = v[i];
        }
        return combined;
    }
    
    // Helper to set state from combined vector
    void set_combined(const Vector& combined) {
        size_t n = u.size();
        for (size_t i = 0; i < n; ++i) {
            u[i] = combined[i];
            v[i] = combined[i + n];
        }
    }
};

} // namespace physics_wave
} // namespace numerical_methods
