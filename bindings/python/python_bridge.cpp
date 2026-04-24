#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "modules/thermodynamics/heat/1d_implicit/model.hpp"
#include "modules/thermodynamics/heat/1d_implicit/state.hpp"
#include "modules/thermodynamics/heat/2d_implicit/model.hpp"
#include "modules/thermodynamics/heat/2d_implicit/state.hpp"
#include "modules/pressure/1d/model.hpp"
#include "modules/pressure/1d/state.hpp"
#include "modules/wave/1d/model.hpp"
#include "modules/wave/1d/state.hpp"
#include "modules/fluids/burgers/model.hpp"
#include "modules/fluids/burgers/state.hpp"
#include "modules/fluids/fluid_dynamics/model.hpp"
#include "modules/fluids/fluid_dynamics/state.hpp"
#include "modules/reservoir/1d/model.hpp"
#include "modules/reservoir/1d/state.hpp"
#include "modules/reservoir/2d/model.hpp"
#include "modules/reservoir/2d/state.hpp"
#include "lib/simulation_engine.hpp"
#include "lib/integrators.hpp"
#include "lib/solvers.hpp"
#include "lib/linearizers.hpp"
#include "lib/engine_infra.hpp"
#include "lib/fem.hpp"

namespace py = pybind11;
using namespace num;
using namespace mod;
using namespace top;

/**
 * @brief Numerical Wrapper for Heat Simulations
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
        model = std::make_shared<heat::Heat1DModel>(cond, storage, TL, TR);
        
        auto discretizer = std::make_shared<heat::Heat1DDiscretizer>();
        auto integrator = std::make_shared<ImplicitEulerIntegrator>();
        auto solver = std::make_shared<LinearTridiagonalSolver>();
        auto linearizer = std::make_shared<NewtonRaphson>(1e-6, 1, false);
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
 * @brief Numerical Wrapper for 2D Heat Simulations
 */
class Heat2DSimulationWrapper {
private:
    std::shared_ptr<heat::Heat2DModel> model;
    std::shared_ptr<heat::Heat2DImplicitState> state;
    std::unique_ptr<SimulationEngine> engine;

public:
    Heat2DSimulationWrapper(int nx, int ny, double Lx, double Ly, double alpha) {
        auto spatial = std::make_shared<Spatial2D>(nx, ny, Lx, Ly);
        state = std::make_shared<heat::Heat2DImplicitState>(spatial, 0.0);
        double dx = Lx / (nx - 1);
        double dy = Ly / (ny - 1);
        auto cond = num::discretization::heat_cond_2d(nx, ny, dx, dy, alpha, 1.0);
        Vector storage(nx * ny, 1.0);
        model = std::make_shared<heat::Heat2DModel>(cond, storage, 0.0, 0.0, 0.0, 0.0);
        
        auto discretizer = std::make_shared<heat::Heat2DDiscretizer>();
        auto integrator = std::make_shared<ImplicitEulerIntegrator>();
        auto solver = std::make_shared<BiCGSTABSolver>();
        auto linearizer = std::make_shared<NewtonRaphson>(1e-6, 1, false);
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
 * @brief Numerical Wrapper for Pressure Simulations
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
        model = std::make_shared<pressure::Pressure1DModel>(cond, storage, PL, PR);
        
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

/**
 * @brief Numerical Wrapper for Wave Simulations
 */
class WaveSimulationWrapper {
private:
    std::shared_ptr<wave::Wave1DModel> model;
    std::shared_ptr<wave::Wave1DState> state;
    std::unique_ptr<SimulationEngine> engine;

public:
    WaveSimulationWrapper(std::shared_ptr<num::discretization::Conductance1D> cond, const Vector& storage) {
        int nx = (int)storage.size();
        auto spatial = std::make_shared<Spatial1D>(nx, 1.0);
        state = std::make_shared<wave::Wave1DState>(spatial, 0.0);
        model = std::make_shared<wave::Wave1DModel>(cond, storage);
        
        auto discretizer = std::make_shared<wave::Wave1DDiscretizer>();
        auto integrator = std::make_shared<ImplicitEulerIntegrator>();
        auto solver = std::make_shared<LinearTridiagonalSolver>();
        auto linearizer = std::make_shared<NewtonRaphson>(1e-6, 1, false);
        auto pm = std::make_shared<SerialParallelManager>();

        engine = std::make_unique<SimulationEngine>(spatial, model, discretizer, integrator, linearizer, solver, pm);
    }

    void set_initial_condition(const std::vector<double>& u, const std::vector<double>& v) {
        state->u = u;
        state->v = v;
    }

    void step(double dt) {
        engine->simulate_step(dt, *state);
    }

    std::vector<double> get_values() const {
        return state->to_vector();
    }
};

/**
 * @brief Numerical Wrapper for Burgers Simulations
 */
class BurgersSimulationWrapper {
private:
    std::shared_ptr<burgers::BurgersModel> model;
    std::shared_ptr<burgers::BurgersState> state;
    std::unique_ptr<SimulationEngine> engine;

public:
    BurgersSimulationWrapper(double nu, double dx, const std::vector<double>& ic) {
        int nx = (int)ic.size();
        auto spatial = std::make_shared<Spatial1D>(nx, dx);
        state = std::make_shared<burgers::BurgersState>(spatial, 0.0);
        state->u = ic;
        model = std::make_shared<burgers::BurgersModel>(nu, dx);
        
        auto discretizer = std::make_shared<burgers::BurgersDiscretizer>();
        auto integrator = std::make_shared<ImplicitEulerIntegrator>();
        auto solver = std::make_shared<LinearTridiagonalSolver>();
        auto linearizer = std::make_shared<NewtonRaphson>(1e-6, 20, true);
        auto pm = std::make_shared<SerialParallelManager>();

        engine = std::make_unique<SimulationEngine>(spatial, model, discretizer, integrator, linearizer, solver, pm);
    }

    void step(double dt) {
        engine->simulate_step(dt, *state);
    }

    std::vector<double> get_values() const {
        return state->to_vector();
    }
};

/**
 * @brief Numerical Wrapper for Stokes Simulations
 */
class StokesSimulationWrapper {
private:
    std::shared_ptr<mod::fluid::FluidModel> model;
    std::shared_ptr<mod::fluid::FluidState> state;
    std::unique_ptr<SimulationEngine> engine;

public:
    StokesSimulationWrapper(int nx, int ny, double Lx, double Ly, double mu) {
        auto mesh = std::make_shared<num::fem::Mesh>();
        mesh->generate_quad_mesh(nx, ny, Lx, Ly);
        state = std::make_shared<mod::fluid::FluidState>(mesh);
        model = std::make_shared<mod::fluid::FluidModel>(mesh, mu, 1.0);
        
        auto discretizer = std::make_shared<mod::fluid::FluidDiscretizer>();
        auto integrator = std::make_shared<ImplicitEulerIntegrator>();
        auto solver = std::make_shared<BiCGSTABSolver>();
        auto linearizer = std::make_shared<NewtonRaphson>(1e-6, 1, false);
        auto pm = std::make_shared<SerialParallelManager>();

        engine = std::make_unique<SimulationEngine>(nullptr, model, discretizer, integrator, linearizer, solver, pm);
    }

    void set_boundary_condition(int node, double u, double v) {
        model->set_velocity_bc(node, u, v);
    }

    void solve() {
        engine->simulate_step(1.0, *state);
    }

    std::vector<double> get_u() const { return state->u; }
    std::vector<double> get_v() const { return state->v; }
    std::vector<double> get_p() const { return state->p; }
};

/**
 * @brief Numerical Wrapper for Reservoir Simulations
 */
class ReservoirSimulationWrapper {
private:
    std::shared_ptr<reservoir::Reservoir1DModel> model;
    std::shared_ptr<reservoir::Reservoir1DState> state;
    std::unique_ptr<SimulationEngine> engine;

public:
    ReservoirSimulationWrapper(std::shared_ptr<num::discretization::Conductance1D> cond, const Vector& storage) {
        int nx = (int)storage.size();
        auto spatial = std::make_shared<Spatial1D>(nx, 1.0);
        state = std::make_shared<reservoir::Reservoir1DState>(spatial, 0.0);
        std::vector<std::shared_ptr<ISourceSink>> wells;
        model = std::make_shared<reservoir::Reservoir1DModel>(cond, storage, wells);
        
        auto discretizer = std::make_shared<reservoir::Reservoir1DDiscretizer>();
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

/**
 * @brief Numerical Wrapper for 2D Reservoir Simulations
 */
class Reservoir2DSimulationWrapper {
private:
    std::shared_ptr<reservoir::Reservoir2DModel> model;
    std::shared_ptr<reservoir::Reservoir2DState> state;
    std::unique_ptr<SimulationEngine> engine;

public:
    Reservoir2DSimulationWrapper(int nx, int ny, double Lx, double Ly, double eta) {
        auto spatial = std::make_shared<Spatial2D>(nx, ny, Lx, Ly);
        state = std::make_shared<reservoir::Reservoir2DState>(spatial, 0.0);
        double dx = Lx / (nx - 1);
        double dy = Ly / (ny - 1);
        auto cond = num::discretization::reservoir_cond_2d(nx, ny, dx, dy, 100.0, 1.0, 1.0, 100.0); // Simplified
        Vector storage(nx * ny, 0.001);

        std::vector<std::shared_ptr<ISourceSink>> wells;
        model = std::make_shared<reservoir::Reservoir2DModel>(cond, storage, wells);
        
        auto discretizer = std::make_shared<reservoir::Reservoir2DDiscretizer>();
        auto integrator = std::make_shared<ImplicitEulerIntegrator>();
        auto solver = std::make_shared<BiCGSTABSolver>();
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

    py::class_<Heat2DSimulationWrapper>(m, "Heat2D")
        .def(py::init<int, int, double, double, double>())
        .def("set_initial_condition", &Heat2DSimulationWrapper::set_initial_condition)
        .def("step", &Heat2DSimulationWrapper::step)
        .def("get_values", &Heat2DSimulationWrapper::get_values);

    py::class_<PressureSimulationWrapper>(m, "Pressure1D")
        .def(py::init([](const std::vector<double>& T_cond, const std::vector<double>& storage, double PL, double PR) {
            auto cond = std::make_shared<num::discretization::Conductance1D>(T_cond.size() + 1);
            cond->T = T_cond;
            return new PressureSimulationWrapper(cond, storage, PL, PR);
        }))
        .def("set_initial_condition", &PressureSimulationWrapper::set_initial_condition)
        .def("step", &PressureSimulationWrapper::step)
        .def("get_values", &PressureSimulationWrapper::get_values);

    py::class_<WaveSimulationWrapper>(m, "Wave1D")
        .def(py::init([](const std::vector<double>& T_cond, const std::vector<double>& storage) {
            auto cond = std::make_shared<num::discretization::Conductance1D>(T_cond.size() + 1);
            cond->T = T_cond;
            return new WaveSimulationWrapper(cond, storage);
        }))
        .def("set_initial_condition", &WaveSimulationWrapper::set_initial_condition)
        .def("step", &WaveSimulationWrapper::step)
        .def("get_values", &WaveSimulationWrapper::get_values);

    py::class_<BurgersSimulationWrapper>(m, "Burgers1D")
        .def(py::init<double, double, const std::vector<double>&>())
        .def("step", &BurgersSimulationWrapper::step)
        .def("get_values", &BurgersSimulationWrapper::get_values);

    py::class_<StokesSimulationWrapper>(m, "Stokes2D")
        .def(py::init<int, int, double, double, double>())
        .def("set_boundary_condition", &StokesSimulationWrapper::set_boundary_condition)
        .def("solve", &StokesSimulationWrapper::solve)
        .def("get_u", &StokesSimulationWrapper::get_u)
        .def("get_v", &StokesSimulationWrapper::get_v)
        .def("get_p", &StokesSimulationWrapper::get_p);

    py::class_<ReservoirSimulationWrapper>(m, "Reservoir1D")
        .def(py::init([](const std::vector<double>& T_cond, const std::vector<double>& storage) {
            auto cond = std::make_shared<num::discretization::Conductance1D>(T_cond.size() + 1);
            cond->T = T_cond;
            return new ReservoirSimulationWrapper(cond, storage);
        }))
        .def("set_initial_condition", &ReservoirSimulationWrapper::set_initial_condition)
        .def("step", &ReservoirSimulationWrapper::step)
        .def("get_values", &ReservoirSimulationWrapper::get_values);

    py::class_<Reservoir2DSimulationWrapper>(m, "Reservoir2D")
        .def(py::init<int, int, double, double, double>())
        .def("set_initial_condition", &Reservoir2DSimulationWrapper::set_initial_condition)
        .def("step", &Reservoir2DSimulationWrapper::step)
        .def("get_values", &Reservoir2DSimulationWrapper::get_values);
}
