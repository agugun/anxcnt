# 2D Single-Phase Reservoir Simulation

## Scenario
This case models a 2D horizontal reservoir drainage. A single producer is perforated at the center of a square drainage area. The simulation visualizes the "cone of depression" as pressure spreads radially from the producer block until it hits the reservoir boundaries.

## Mathematics
The model solves the 2D pressure diffusivity equation:
$$\frac{\partial p}{\partial t} = \eta \left( \frac{\partial^2 p}{\partial x^2} + \frac{\partial^2 p}{\partial y^2} \right) - \text{Wells}$$

The spatial discretization uses a standard **5-point central difference stencil**.

## Logic
- **Solver**: Uses the **Conjugate Gradient (CG)** iterative solver.
- **Matrix-Free**: Instead of building a large 2D sparse matrix, the model implements an analytical `apply_jacobian` (Matrix-vector product). This minimizes memory overhead and allows scaling to larger grids.
- **Boundaries**: **No-flow** boundaries on all four sides.
- **Visuals**: Multi-field VTK output (`.vti`) allows viewing `Pressure` and `WellLocation` simultaneously in ParaView.
