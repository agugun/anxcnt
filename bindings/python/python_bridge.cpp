#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "modules/thermodynamics/heat/1d_implicit/model.hpp"
#include "modules/thermodynamics/heat/1d_implicit/state.hpp"
#include "modules/pressure/1d/model.hpp"
#include "modules/pressure/1d/state.hpp"
#include "lib/simulation_engine.hpp"
#include "lib/integrators.hpp"
#include "lib/solvers.hpp"
#include "lib/linearizers.hpp"
#include "lib/engine_infra.hpp"

namespace py = pybind11;
using namespace num;
using namespace mod;
using namespace top;

/**
 * @brief Numerical Wrapper for Heat Simulations (Calculus in Python)
 */
class HeatSimulationWrapper {
private:
    std::shared_ptr<heat::Heat1DModel> model;
    std::shared_ptr<heat::Heat1DImplicitState> state;
    std::unique_ptr<SimulationEngine> engine;

public:
    HeatSimulationWrapper(std::shared_ptr<num::discretization::Conductance1D> cond, const Vector& storage, double TL, double TR) {
        int nx = (int)storage.size();
        auto spatial = std::make_shared<Spatial1D>(nx, 1.0);
        state = std::make_shared<heat::Heat1DImplicitState>(spatial, 0.0);
        
        Vector storage_coeffs = storage;
        model = std::make_shared<heat::Heat1DModel>(cond, storage_coeffs, TL, TR);
        
        auto discretizer = std::make_shared<heat::Heat1DDiscretizer>();
        auto integrator = std::make_shared<ImplicitEulerIntegrator>();
        auto solver = std::make_shared<LinearTridiagonalSolver>();
        auto linearizer = std::make_shared<NewtonRaphson>(1e-6, 1, false); // Linear system
        auto pm = std::make_shared<SerialParallelManager>();

        engine = std::make_unique<SimulationEngine>(spatial, model, discretizer, integrator, linearizer, solver, pm);
    }

    void set_initial_condition(const std::vector<double>& ic) {
        state->temperatures = ic;
    }

    void step(double dt) {
        engine->simulate_step(dt, *state);
    }

    std::vector<double> get_values() const {
        return state->to_vector();
    }
};

/**
 * @brief Numerical Wrapper for Pressure Simulations (Calculus in Python)
 */
class PressureSimulationWrapper {
private:
    std::shared_ptr<pressure::Pressure1DModel> model;
    std::shared_ptr<pressure::Pressure1DState> state;
    std::unique_ptr<SimulationEngine> engine;

public:
    PressureSimulationWrapper(std::shared_ptr<num::discretization::Conductance1D> cond, const Vector& storage, double PL, double PR) {
        int nx = (int)storage.size();
        auto spatial = std::make_shared<Spatial1D>(nx, 1.0);
        state = std::make_shared<pressure::Pressure1DState>(spatial, 0.0);
        
        Vector storage_coeffs = storage;
        model = std::make_shared<pressure::Pressure1DModel>(cond, storage_coeffs, PL, PR);
        
        auto discretizer = std::make_shared<pressure::Pressure1DDiscretizer>();
        auto integrator = std::make_shared<ImplicitEulerIntegrator>();
        auto solver = std::make_shared<LinearTridiagonalSolver>();
        auto linearizer = std::make_shared<NewtonRaphson>(1e-6, 1, false);
        auto pm = std::make_shared<SerialParallelManager>();

        engine = std::make_unique<SimulationEngine>(spatial, model, discretizer, integrator, linearizer, solver, pm);
    }

    void set_initial_condition(const std::vector<double>& ic) {
        state->pressures = ic;
    }

    void step(double dt) {
        engine->simulate_step(dt, *state);
    }

    std::vector<double> get_values() const {
        return state->to_vector();
    }
};

PYBIND11_MODULE(axcnt_cpp, m) {
    m.doc() = "NumPhys Core Python Bridge - Numerical Execution Engine";

    py::class_<num::discretization::Conductance1D, std::shared_ptr<num::discretization::Conductance1D>>(m, "Conductance1D")
        .def(py::init<size_t>())
        .def_readwrite("T", &num::discretization::Conductance1D::T);

    py::class_<HeatSimulationWrapper>(m, "Heat1D")
        .def(py::init([](const std::vector<double>& T_cond, const std::vector<double>& storage, double TL, double TR) {
            auto cond = std::make_shared<num::discretization::Conductance1D>(T_cond.size() + 1);
            cond->T = T_cond;
            return new HeatSimulationWrapper(cond, storage, TL, TR);
        }))
        .def("set_initial_condition", &HeatSimulationWrapper::set_initial_condition)
        .def("step", &HeatSimulationWrapper::step)
        .def("get_values", &HeatSimulationWrapper::get_values);

    py::class_<PressureSimulationWrapper>(m, "Pressure1D")
        .def(py::init([](const std::vector<double>& T_cond, const std::vector<double>& storage, double PL, double PR) {
            auto cond = std::make_shared<num::discretization::Conductance1D>(T_cond.size() + 1);
            cond->T = T_cond;
            return new PressureSimulationWrapper(cond, storage, PL, PR);
        }))
        .def("set_initial_condition", &PressureSimulationWrapper::set_initial_condition)
        .def("step", &PressureSimulationWrapper::step)
        .def("get_values", &PressureSimulationWrapper::get_values);
}
