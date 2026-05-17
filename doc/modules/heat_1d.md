# 1D Heat Conduction Module (Implicit)

Simulates transient heat transfer in a 1D rod using an Implicit Euler scheme.

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

    namespace Heat_1D_Module {
        class Heat1DImplicitState {
            <<IState>>
            +temperatures: Vector
            +spatial: shared_ptr~Spatial1D~
            +clone() unique_ptr
        }
        class Heat1DModel {
            <<IModel>>
            +T_left: double
            +T_right: double
            +cond: shared_ptr~Conductance1D~
            +storage_coeff: Vector
            +build_capacity(grd: IGrid, st: IState) Vector
        }
        class Heat1DDiscretizer {
            <<IDiscretizer>>
            +build_jacobian(grd: IGrid, mdl: IModel, st: IState, J: SparseMatrix)
            +build_residual(grd: IGrid, mdl: IModel, st: IState, R: Vector)
            +apply_bc(grd: IGrid, mdl: IModel, st: IState, J: SparseMatrix, R: Vector)
        }
        class HeatSimulation {
            <<Simulation>>
            +build(config: ConfigReader)
            +create_initial_state(config: ConfigReader) unique_ptr
        }
    }

    %% Relationships
    Simulation o-- IModel : Aggregation
    Simulation o-- IDiscretizer : Aggregation
    Simulation o-- IState : Aggregation
    
    HeatSimulation ..> Heat1DModel : Dependency (uses)
    HeatSimulation ..> Heat1DImplicitState : Dependency (uses)
    HeatSimulation ..> Heat1DDiscretizer : Dependency (uses)
```
