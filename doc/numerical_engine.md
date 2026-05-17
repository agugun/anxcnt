# Numerical Engine Architecture (`num`)

This document details the numerical strategies and solver implementations that power the AXSCNT simulation engine. The core philosophy is **Numerical Strategy Injection**, where solvers and integrators are treated as swappable components.

## Class Diagram

```mermaid
classDiagram
    direction TB

    namespace Framework_Interfaces {
        class ITimeIntegrator {
            <<interface>>
            +apply_temporal(grd, mdl, J, R, st_new, st_old, dt)*
            +compute_dt(st, t)*
        }
        class ILinearizer {
            <<interface>>
            +resolve(st_n, dt, grd, mdl, disc, timer, solver, parallel)*
        }
        class ISolver {
            <<interface>>
            +solve(A, b)*
        }
    }

    namespace Integrators {
        class ImplicitEulerIntegrator {
            <<ITimeIntegrator>>
            +apply_temporal()
        }
        class ExplicitEulerIntegrator {
            <<ITimeIntegrator>>
            +apply_temporal()
        }
        class RungeKutta4Integrator {
            <<ITimeIntegrator>>
            +apply_temporal()
        }
    }

    namespace Linearizers {
        class NewtonRaphson {
            <<ILinearizer>>
            -tolerance: double
            -max_iterations: int
            +resolve()
        }
        class PicardIteration {
            <<ILinearizer>>
            +resolve()
        }
    }

    namespace Solvers {
        class LUSolver {
            <<ISolver>>
            +solve(A, b)
        }
        class BiCGSTABSolver {
            <<ISolver>>
            +solve(A, b)
        }
        class LinearTridiagonalSolver {
            <<ISolver>>
            +solve(A, b)
        }
    }

    namespace Utility_Boundary {
        class IParallelManager {
            <<utl>>
            +sync_ghost_cells(st)
            +get_global_norm(r)
        }
        class SerialParallelManager {
            <<utl::IParallelManager>>
        }
    }

    %% Relationships (Contract Fulfillment)
    ITimeIntegrator <|.. ImplicitEulerIntegrator
    ITimeIntegrator <|.. ExplicitEulerIntegrator

    ILinearizer <|.. NewtonRaphson

    ISolver <|.. LUSolver
    ISolver <|.. BiCGSTABSolver
    ISolver <|.. LinearTridiagonalSolver

    IParallelManager <|.. SerialParallelManager
```

`IParallelManager` and `SerialParallelManager` are utility-layer hooks in `namespace utl`, but the numerical linearizers receive them so convergence checks and future distributed halo synchronization can use the same orchestration path.

## Numerical Strategies

### ⏱️ Time Integration
The engine supports both **Implicit** and **Explicit** schemes.
*   **Implicit Euler**: Preferred for stiff equations (like Heat and Pressure) as it is unconditionally stable.
*   **Explicit schemes**: Available for wave propagation or non-stiff dynamics.

### 📉 Linearization (Newton-Raphson)
For non-linear physics, the `NewtonRaphson` class implements the iterative process:
1.  **Assembly**: Calls `IDiscretizer` to build the Jacobian ($J$) and Residual ($R$).
2.  **Temporal Correction**: Calls `ITimeIntegrator` to inject mass-matrix terms.
3.  **Solve**: Calls `ISolver` to find the update $\delta = -J^{-1}R$.
4.  **Update**: Calls `IState` to apply the correction.

### 🧮 Linear Solvers
*   **Tridiagonal**: Optimized $O(N)$ solver for 1D structured problems.
*   **BiCGSTAB**: Iterative solver for large, sparse 2D/3D systems.
*   **LU**: Direct solver for small, dense, or specific sparse matrices.
