#pragma once
#include <memory>
#include <vector>
#include "../../lib/model/state/field.hpp"
#include "../../lib/simulation/newton_simulator.hpp"
#include "../../lib/numerical_methods/solver/linear_solver.hpp"
#include "model.hpp"

namespace numerical_methods::physics_pressure {

class PressureSimulation : public simulation::NewtonSimulator {
public:
    PressureSimulation(int nx, double dx, double k, double phi, double mu, double ct)
        : state(nx, dx) {
        problem = std::make_shared<model::Pressure1DModel>(k, phi, mu, ct, 3000.0, 3000.0);
        this->set_solver(std::make_shared<numerical_methods::solver::LinearSolver>());
    }

    void step(double dt) {
        simulation::NewtonSimulator::step(dt, *problem, state);
        problem->applyBoundaryConditions(state);
    }

    void set_initial_condition(const std::vector<double>& ic) {
        state.set_initial_conditions(ic);
    }

    void set_boundary_conditions(double left, double right) {
        problem->set_boundary_conditions({left, right});
        problem->applyBoundaryConditions(state);
    }

    const std::vector<double>& get_values() const {
        return state.get_data();
    }

    const model::state::Field& getState() const {
        return state;
    }

private:
    model::state::Field state;
    std::shared_ptr<model::Pressure1DModel> problem;
};

} // namespace numerical_methods::physics_pressure
