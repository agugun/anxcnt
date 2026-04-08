#pragma once
#include "../../model/state/field.hpp"
#include "../../model/state/vector_field.hpp"

namespace numerical_methods::operator_unit {

/**
 * @brief Gradient operator.
 * 
 * Maps Cell-centered values to Face-centered values.
 */
class Gradient {
public:
    /**
     * @brief Computes the gradient of a scalar field.
     * 
     * @param input Scalar field at cell centers.
     * @param output Vector field at cell faces (must have size input.nx + 1).
     */
    void apply(const model::state::Field& input, model::state::VectorField& output) const {
        int nx = input.nx;
        double dx = input.dx;

        // Internal faces
        for (int i = 1; i < nx; ++i) {
            output[i] = (input[i] - input[i-1]) / dx;
        }

        // Boundary faces (assuming zero gradient placeholders if not specified)
        output[0] = 0.0;
        output[nx] = 0.0;
    }
};

} // namespace numerical_methods::operator_unit
