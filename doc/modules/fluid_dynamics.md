# 2D Fluid Dynamics Module (FEM)

This module implements an incompressible Navier-Stokes solver using the Finite Element Method (FEM) on unstructured triangular meshes.

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

    namespace Fluid_Module {
        class FluidState {
            <<IState>>
            +u: Vector
            +v: Vector
            +p: Vector
            +mesh: shared_ptr~Mesh~
            +update(delta: Vector)
            +clone() unique_ptr
        }
        class FluidModel {
            <<IModel>>
            +mu: double
            +rho: double
            +mesh: shared_ptr~Mesh~
            +velocity_bc_nodes: vector~int~
            +set_velocity_bc(node_idx, u_val, v_val)
            +build_capacity(grd: IGrid, st: IState) Vector
        }
        class FluidDiscretizer {
            <<IDiscretizer>>
            +build_jacobian(grd: IGrid, mdl: IModel, st: IState, J: SparseMatrix)
            +build_residual(grd: IGrid, mdl: IModel, st: IState, R: Vector)
            +apply_bc(grd: IGrid, mdl: IModel, st: IState, J: SparseMatrix, R: Vector)
        }
        class FluidSimulation {
            <<Simulation>>
            +build(config: ConfigReader)
            +create_initial_state(config: ConfigReader) unique_ptr
        }
    }

    %% Relationships
    Simulation o-- IModel : Aggregation
    Simulation o-- IDiscretizer : Aggregation
    Simulation o-- IState : Aggregation
    
    FluidSimulation ..> FluidModel : Dependency (uses)
    FluidSimulation ..> FluidState : Dependency (uses)
    FluidSimulation ..> FluidDiscretizer : Dependency (uses)
```
