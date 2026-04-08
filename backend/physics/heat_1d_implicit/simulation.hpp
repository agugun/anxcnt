#pragma once
#include <memory>
#include <vector>
#include "../../lib/model/state/field.hpp"
#include "../../lib/simulation/transient_simulator.hpp"
#include "../../lib/numerical_methods/solver/linear_solver.hpp"
#include "model.hpp"

namespace numerical_methods::physics_heat {

class HeatSimulation : public simulation::TransientSimulator {
public:
    HeatSimulation(int nx, double dx, double alpha)
        : state(nx, dx) {
        problem = std::make_shared<model::Heat1DModel>(alpha, 0.0, 0.0);
        this->set_solver(std::make_shared<numerical_methods::solver::LinearSolver>());
    }

    void step(double dt) {
        // Utilize the base transient simulator step
        simulation::TransientSimulator::step(dt, *problem, state);
        
        // Ensure boundary conditions are applied to the state explicitly if needed
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
    std::shared_ptr<model::Heat1DModel> problem;
};

} // namespace numerical_methods::physics_heat
