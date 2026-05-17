# Harmonic Oscillator Module

The Harmonic Oscillator serves as the "Hello World" of AXSCNT, demonstrating the 0D ODE integration within the framework.

## Class Diagram

```mermaid
classDiagram
    direction TB

    namespace Framework_Contracts {
        class IState {
            <<interface>>
            +update(delta: Vector)*
            +to_vector() Vector*
            +clone() unique_ptr~IState~*
        }
        class IModel {
            <<interface>>
            +get_tolerance() double*
            +build_capacity(grd: IGrid, st: IState) Vector*
        }
        class IDiscretizer {
            <<interface>>
            +build_jacobian(grd: IGrid, mdl: IModel, st: IState, J: SparseMatrix)*
            +build_residual(grd: IGrid, mdl: IModel, st: IState, R: Vector)*
            +apply_bc(grd: IGrid, mdl: IModel, st: IState, J: SparseMatrix, R: Vector)*
        }
        class Simulation {
            <<abstract>>
            #grd: shared_ptr~IGrid~
            #mdl: shared_ptr~IModel~
            #timer: shared_ptr~ITimeIntegrator~
            #linearizer: shared_ptr~ILinearizer~
            #solver: shared_ptr~ISolver~
            +build(config: ConfigReader)*
            +run(t_max: double, dt: double, st_init: unique_ptr) unique_ptr
        }
    }

    namespace Oscillator_Module {
        class OscillatorState {
            <<IState>>
            +x: Position_m
            +v: Velocity_ms
            +update(delta: Vector)
            +clone() unique_ptr
        }
        class OscillatorModel {
            <<IModel>>
            +m: Mass_kg
            +c: Damping_Ns_m
            +k: Stiffness_N_m
            +build_capacity(grd: IGrid, st: IState) Vector
        }
        class OscillatorDiscretizer {
            <<IDiscretizer>>
            +build_jacobian(grd: IGrid, mdl: IModel, st: IState, J: SparseMatrix)
            +build_residual(grd: IGrid, mdl: IModel, st: IState, R: Vector)
            +apply_bc(grd: IGrid, mdl: IModel, st: IState, J: SparseMatrix, R: Vector)
        }
        class OscillatorSimulation {
            <<Simulation>>
            +build(config: ConfigReader)
            +create_initial_state(config: ConfigReader) unique_ptr
        }
    }

    %% Relationships
    Simulation o-- IModel : Aggregation
    Simulation o-- IDiscretizer : Aggregation
    Simulation o-- IState : Aggregation
    
    OscillatorSimulation ..> OscillatorModel : Dependency (uses)
    OscillatorSimulation ..> OscillatorState : Dependency (uses)
    OscillatorSimulation ..> OscillatorDiscretizer : Dependency (uses)
```

## Physics Equations

- $dx/dt = v$
- $m \cdot dv/dt + c \cdot v + k \cdot x = 0$
