#pragma once
#include "../../model/state/field.hpp"
#include "../../model/state/vector_field.hpp"

namespace numerical_methods::operator_unit {

/**
 * @brief Divergence operator.
 * 
 * Maps Face-centered values (fluxes) to Cell-centered values.
 */
class Divergence {
public:
    /**
     * @brief Computes the divergence of a vector field.
     * 
     * @param input Vector field at cell faces.
     * @param output Scalar field at cell centers (must have size input.nfaces - 1).
     */
    void apply(const model::state::VectorField& input, model::state::Field& output) const {
        int nx = output.nx;
        double dx = output.dx;

        for (int i = 0; i < nx; ++i) {
            output[i] = (input[i+1] - input[i]) / dx;
        }
    }
};

} // namespace numerical_methods::operator_unit
