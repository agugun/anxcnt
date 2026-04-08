#pragma once
#include "i_operator.hpp"
#include "gradient.hpp"
#include "divergence.hpp"

namespace numerical_methods::operator_unit {

/**
 * @brief Laplacian operator (div grad).
 * 
 * Implements the second-order central difference scheme as a composition
 * of gradient (Cell->Face) and divergence (Face->Cell).
 */
class Laplacian : public IOperator {
public:
    void apply(const model::state::Field& input, std::vector<double>& output) const override {
        int nx = input.nx;
        model::state::Field output_field(nx, input.dx); // Temporary field to use with Divergence
        model::state::VectorField grad(nx + 1, input.dx);
        
        Gradient gradient;
        Divergence divergence;

        gradient.apply(input, grad);
        divergence.apply(grad, output_field);

        // Copy all points of the result to the output vector
        for (int i = 0; i < nx; ++i) {
            output[i] = output_field[i];
        }
    }

    void assembleWithDx(std::vector<double>& a, 
                  std::vector<double>& b, 
                  std::vector<double>& c, 
                  std::vector<double>& d,
                  double dx,
                  double scale = 1.0) const override {
        int m = b.size();
        double dx2 = dx * dx;
        double coeff = scale / dx2;

        for (int i = 0; i < m; ++i) {
            b[i] += -2.0 * coeff;
            if (i > 0) a[i] += coeff;
            if (i < m - 1) c[i] += coeff;
        }
    }
};

} // namespace numerical_methods::operator_unit
