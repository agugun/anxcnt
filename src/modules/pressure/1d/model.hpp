#pragma once
#include "lib/interfaces.hpp"
#include "state.hpp"
#include "lib/discretization.hpp"
 
namespace mod {
using namespace top;
namespace pressure {

/**
 * @brief 1D Pressure Physical Model (Properties).
 */
class Pressure1DModel : public IModel {
public:
    double p_left, p_right;
    std::shared_ptr<num::discretization::Conductance1D> cond;
    Vector storage_coeff;

    Pressure1DModel(std::shared_ptr<num::discretization::Conductance1D> c, const Vector& storage, double pl, double pr)
        : p_left(pl), p_right(pr), cond(c), storage_coeff(storage) {}

    double get_tolerance() const override { return 1e-6; }

    Vector get_accumulation_weights(const IGrid& grid, const IState& state) const override {
        return storage_coeff;
    }
};


} // namespace pressure
} // namespace mod
