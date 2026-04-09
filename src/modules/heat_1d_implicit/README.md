# Heat Equation 1D (Implicit)

## Physics Background
This module simulates 1D heat conduction using an implicit time-integration scheme (Implicit Euler). Implicit methods are unconditionally stable, allowing for much larger time steps than explicit methods.

## Governing Equation
$$\frac{\partial T}{\partial t} = \alpha \frac{\partial^2 T}{\partial x^2}$$

## Numerical Implementation
- **Time Integration**: Implicit Euler:
  $$\frac{T^{n+1} - T^n}{\Delta t} = \alpha \frac{\partial^2 T^{n+1}}{\partial x^2}$$
- **Linear System**: The discretization leads to a tridiagonal matrix system:
  $$(I - \Delta t \cdot \alpha \Delta) T^{n+1} = T^n$$
  $$-r T_{i-1}^{n+1} + (1 + 2r) T_i^{n+1} - r T_{i+1}^{n+1} = T_i^n$$
  where $r = \frac{\alpha \Delta t}{\Delta x^2}$.
- **Solver**: Uses a specialized **Linear Tridiagonal Solver** (Thomas Algorithm) for $O(N)$ efficiency.
