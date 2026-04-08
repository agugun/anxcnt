#pragma once
#include <vector>

namespace model {

/**
 * @brief Data structure representing the assembled algebraic system.
 * 
 * For now, mainly used passing tridiagonal matrices (a, b, c) and 
 * Right-Hand Side vectors (d) between IModel and ISolver.
 */
struct SystemMatrix {
    std::vector<double> a; // Lower diagonal
    std::vector<double> b; // Main diagonal
    std::vector<double> c; // Upper diagonal
    std::vector<double> d; // Right-hand side

    void resize(int size) {
        a.resize(size, 0.0);
        b.resize(size, 0.0);
        c.resize(size, 0.0);
        d.resize(size, 0.0);
    }
};

} // namespace model
