#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "physics/heat_1d_implicit/simulation.hpp"
#include "physics/pressure_diffusivity_1d/simulation.hpp"

namespace py = pybind11;

PYBIND11_MODULE(cnt, m) {
    m.doc() = "NumPhys Backend Python Bridge";

    // Bind the new Case-based simulation orchestrator
    py::class_<numerical_methods::physics_heat::HeatSimulation>(m, "Heat1DImplicit")
        .def(py::init<int, double, double>(), 
             py::arg("nx"), py::arg("dx"), py::arg("alpha"))
        .def("set_initial_condition", &numerical_methods::physics_heat::HeatSimulation::set_initial_condition,
             py::arg("ic"))
        .def("set_boundary_conditions", &numerical_methods::physics_heat::HeatSimulation::set_boundary_conditions,
             py::arg("left"), py::arg("right"))
        .def("step", &numerical_methods::physics_heat::HeatSimulation::step,
             py::arg("dt"))
        .def("get_values", [](const numerical_methods::physics_heat::HeatSimulation& self) {
            return self.get_values();
        });

    // Bind the new Pressure Diffusivity Case (Newton-Raphson)
    py::class_<numerical_methods::physics_pressure::PressureSimulation>(m, "Pressure1DImplicit")
        .def(py::init<int, double, double, double, double, double>(),
             py::arg("nx"), py::arg("dx"), py::arg("k"), py::arg("phi"), py::arg("mu"), py::arg("ct"))
        .def("set_initial_condition", &numerical_methods::physics_pressure::PressureSimulation::set_initial_condition,
             py::arg("ic"))
        .def("set_boundary_conditions", &numerical_methods::physics_pressure::PressureSimulation::set_boundary_conditions,
             py::arg("left"), py::arg("right"))
        .def("step", &numerical_methods::physics_pressure::PressureSimulation::step,
             py::arg("dt"))
        .def("get_values", [](const numerical_methods::physics_pressure::PressureSimulation& self) {
            return self.get_values();
        });
}
