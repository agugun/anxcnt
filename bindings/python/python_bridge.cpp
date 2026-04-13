#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "modules/heat/1d_implicit/model.hpp"
#include "modules/heat/1d_implicit/state.hpp"
#include "modules/pressure/1d/model.hpp"
#include "modules/pressure/1d/state.hpp"
#include "lib/integrators.hpp"
#include "lib/spatial.hpp"

namespace py = pybind11;
using namespace num;
using namespace mod;
using namespace top;


/**
 * @brief Helper class to maintain the expected Python API for Heat simulations
 */
class HeatSimulationWrapper {
private:
    std::shared_ptr<mod::heat::Heat1DModel> model;
    std::shared_ptr<mod::heat::Heat1DImplicitState> state;
    std::shared_ptr<LinearTridiagonalSolver> solver;
    std::shared_ptr<ImplicitEulerIntegrator> integrator;
    std::unique_ptr<StandardSimulator> simulator;

public:
    HeatSimulationWrapper(int nx, double dx, double alpha) {
        Spatial1D spatial(nx, dx);
        state = std::make_shared<mod::heat::Heat1DImplicitState>(spatial, 0.0);
        model = std::make_shared<mod::heat::Heat1DModel>(alpha, 0.0, 0.0);
        solver = std::make_shared<LinearTridiagonalSolver>();
        integrator = std::make_shared<ImplicitEulerIntegrator>();
        simulator = std::make_unique<StandardSimulator>(model, state, solver, integrator);
    }

    void set_initial_condition(const std::vector<double>& ic) {
        state->temperatures = ic;
    }

    void set_boundary_conditions(double left, double right) {
        model->set_bcs(left, right);
    }

    void step(double dt) {
        // Step once
        integrator->step(*model, *state, solver.get(), dt);
    }

    std::vector<double> get_values() const {
        return state->temperatures;
    }
};

/**
 * @brief Helper class to maintain the expected Python API for Pressure simulations
 */
class PressureSimulationWrapper {
private:
    std::shared_ptr<mod::pressure::Pressure1DModel> model;
    std::shared_ptr<mod::pressure::Pressure1DState> state;
    std::shared_ptr<LinearTridiagonalSolver> solver;
    std::shared_ptr<ImplicitEulerIntegrator> integrator;
    std::unique_ptr<StandardSimulator> simulator;

public:
    PressureSimulationWrapper(int nx, double dx, double k, double phi, double mu, double ct) {
        Spatial1D spatial(nx, dx);
        state = std::make_shared<mod::pressure::Pressure1DState>(spatial, 0.0);
        model = std::make_shared<mod::pressure::Pressure1DModel>(k, phi, mu, ct, 0.0, 0.0);
        solver = std::make_shared<LinearTridiagonalSolver>();
        integrator = std::make_shared<ImplicitEulerIntegrator>();
        simulator = std::make_unique<StandardSimulator>(model, state, solver, integrator);
    }

    void set_initial_condition(const std::vector<double>& ic) {
        state->pressures = ic;
    }

    void set_boundary_conditions(double left, double right) {
        model->set_bcs(left, right);
    }

    void step(double dt) {
        integrator->step(*model, *state, solver.get(), dt);
    }

    std::vector<double> get_values() const {
        return state->pressures;
    }
};

PYBIND11_MODULE(cnt, m) {
    m.doc() = "NumPhys Core Python Bridge";

    py::class_<HeatSimulationWrapper>(m, "Heat1DImplicit")
        .def(py::init<int, double, double>(), 
             py::arg("nx"), py::arg("dx"), py::arg("alpha"))
        .def("set_initial_condition", &HeatSimulationWrapper::set_initial_condition)
        .def("set_boundary_conditions", &HeatSimulationWrapper::set_boundary_conditions)
        .def("step", &HeatSimulationWrapper::step)
        .def("get_values", &HeatSimulationWrapper::get_values);

    py::class_<PressureSimulationWrapper>(m, "Pressure1DImplicit")
        .def(py::init<int, double, double, double, double, double>(),
             py::arg("nx"), py::arg("dx"), py::arg("k"), py::arg("phi"), py::arg("mu"), py::arg("ct"))
        .def("set_initial_condition", &PressureSimulationWrapper::set_initial_condition)
        .def("set_boundary_conditions", &PressureSimulationWrapper::set_boundary_conditions)
        .def("step", &PressureSimulationWrapper::step)
        .def("get_values", &PressureSimulationWrapper::get_values);
}
