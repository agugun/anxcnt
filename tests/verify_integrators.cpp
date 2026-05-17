#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

#include "lib/engine_infra.hpp"
#include "lib/integrators.hpp"
#include "lib/linearizers.hpp"
#include "lib/simulation.hpp"
#include "lib/solvers.hpp"

using namespace mod;
using namespace num;
using namespace std;

class ScalarGrid : public IGrid {
public:
    size_t get_total_cells() const override { return 1; }
};

class ScalarState : public IState {
public:
    Vector value;

    explicit ScalarState(double v) : value{v} {}

    void update(const Vector& delta) override {
        value[0] += delta[0];
    }

    Vector to_vector() const override {
        return value;
    }

    unique_ptr<IState> clone() const override {
        return make_unique<ScalarState>(value[0]);
    }
};

class UnitCapacityModel : public IModel {
public:
    double get_tolerance() const override { return 1e-12; }

    Vector build_capacity(const IGrid& grd, const IState& st) const override {
        return {1.0};
    }
};

class DecayDiscretizer : public IDiscretizer {
public:
    void build_jacobian(const IGrid& grd, const IModel& mdl, const IState& st, SparseMatrix& J) const override {
        if (J.rows != 1 || J.cols != 1) J = SparseMatrix(1, 1);
        J.triplets.clear();
        J.triplets.push_back({0, 0, 1.0});
    }

    void build_residual(const IGrid& grd, const IModel& mdl, const IState& st, Vector& R) const override {
        const auto& scalar = static_cast<const ScalarState&>(st);
        R[0] = scalar.value[0];
    }

    void apply_bc(const IGrid& grd, const IModel& mdl, const IState& st, SparseMatrix& J, Vector& R) const override {}
};

void test_implicit_euler_temporal_assembly() {
    cout << "--- Testing Implicit Euler Temporal Assembly ---" << endl;

    ScalarGrid grid;
    UnitCapacityModel model;
    ScalarState old_state(1.0);
    ScalarState new_state(1.2);
    SparseMatrix J(1, 1);
    Vector R(1, 0.0);

    ImplicitEulerIntegrator integrator;
    integrator.apply_temporal(grid, model, J, R, new_state, old_state, 0.1);
    J.compress();

    auto dense = J.to_dense();
    cout << "Residual contribution: " << R[0] << endl;
    cout << "Jacobian diagonal contribution: " << dense[0][0] << endl;

    assert(abs(R[0] - 2.0) < 1e-12);
    assert(abs(dense[0][0] - 10.0) < 1e-12);
}

void test_implicit_euler_decay_step() {
    cout << "\n--- Testing Implicit Euler Decay Step ---" << endl;

    auto grid = make_shared<ScalarGrid>();
    auto model = make_shared<UnitCapacityModel>();
    auto discretizer = make_shared<DecayDiscretizer>();
    auto integrator = make_shared<ImplicitEulerIntegrator>();
    auto linearizer = make_shared<NewtonRaphson>(1e-12, 8, false);
    auto solver = make_shared<LUSolver>();
    auto parallel = make_shared<utl::SerialParallelManager>();

    SimulationEngine engine(grid, model, discretizer, integrator, linearizer, solver, parallel);
    ScalarState state(1.0);

    auto next = engine.step(0.0, 0.1, state, 1);
    double numerical = next->to_vector()[0];
    double expected = 1.0 / 1.1;

    cout << "Numerical: " << numerical << endl;
    cout << "Expected: " << expected << endl;

    assert(abs(numerical - expected) < 1e-12);
}

void test_explicit_integrator_guards() {
    cout << "\n--- Testing Explicit Integrator Guards ---" << endl;

    ScalarGrid grid;
    UnitCapacityModel model;
    ScalarState old_state(1.0);
    ScalarState new_state(1.0);
    SparseMatrix J(1, 1);
    Vector R(1, 0.0);

    ForwardEulerIntegrator forward;
    RungeKutta4Integrator rk4;

    bool forward_threw = false;
    bool rk4_threw = false;

    try {
        forward.apply_temporal(grid, model, J, R, new_state, old_state, 0.1);
    } catch (const runtime_error&) {
        forward_threw = true;
    }

    try {
        rk4.apply_temporal(grid, model, J, R, new_state, old_state, 0.1);
    } catch (const runtime_error&) {
        rk4_threw = true;
    }

    assert(forward_threw);
    assert(rk4_threw);
}

int main() {
    try {
        test_implicit_euler_temporal_assembly();
        test_implicit_euler_decay_step();
        test_explicit_integrator_guards();
        cout << "\nAll Integrator Tests Passed." << endl;
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}
