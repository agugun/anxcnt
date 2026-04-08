#ifndef FIELD_HPP
#define FIELD_HPP

#include <vector>
#include <string>

#include "i_state.hpp"

namespace model::state {

struct Field : public IState {
    int nx;
    double dx;
    std::vector<double> values;
    std::string name;

    Field(int n, double d, std::string nme = "Field")
        : nx(n), dx(d), values(n, 0.0), name(nme) {}

    // Implement IState interface
    void set_initial_conditions(const std::vector<double>& initial_values) override {
        if (initial_values.size() == static_cast<size_t>(nx)) {
            values = initial_values;
        }
    }

    const std::vector<double>& get_data() const override {
        return values;
    }

    double& operator[](int i) { return values[i]; }
    const double& operator[](int i) const { return values[i]; }

    double get_coordinate(int i) const { return i * dx; }

    int size() const { return nx; }
};

} // namespace model::state

#endif // FIELD_HPP
