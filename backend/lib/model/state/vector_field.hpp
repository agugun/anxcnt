#ifndef VECTOR_FIELD_HPP
#define VECTOR_FIELD_HPP

#include <vector>
#include <string>

namespace model::state {

/**
 * @brief Represents a vector field, typically defined at cell faces.
 */
struct VectorField {
    int nfaces;
    double dx;
    std::vector<double> values;
    std::string name;

    VectorField(int n, double d, std::string nme = "VectorField")
        : nfaces(n), dx(d), values(n, 0.0), name(nme) {}

    double& operator[](int i) { return values[i]; }
    const double& operator[](int i) const { return values[i]; }

    int size() const { return nfaces; }
};

} // namespace model::state

#endif // VECTOR_FIELD_HPP
