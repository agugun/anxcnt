#pragma once
#include "../../model/state/field.hpp"
#include <vector>

namespace numerical_methods::operator_unit {

/**
 * @brief Base interface for mathematical operators.
 * 
 * Note: Namespace is 'operator_unit' because 'operator' is a C++ keyword.
 */
class IOperator {
public:
    virtual ~IOperator() = default;

    /**
     * @brief Applies the operator to a constant field.
     */
    virtual void apply(const model::state::Field& input, std::vector<double>& output) const = 0;

    /**
     * @brief Assembles the tridiagonal coefficients (for implicit schemes).
     * 
     * @param a Lower diagonal
     * @param b Main diagonal
     * @param c Upper diagonal
     * @param d Residual/RHS vector
     * @param dx Grid spacing
     * @param scale Scaling factor (e.g. dt * diffusivity)
     */
    virtual void assembleWithDx(std::vector<double>& a, 
                                std::vector<double>& b, 
                                std::vector<double>& c, 
                                std::vector<double>& d,
                                double dx,
                                double scale = 1.0) const = 0;
};

} // namespace numerical_methods::operator_unit
