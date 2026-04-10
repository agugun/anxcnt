# Heat Equation 1D (Explicit)

## Physics Background
This module simulates 1D heat conduction in a solid rod using an explicit time-integration scheme (Forward Euler). It models the spatial and temporal evolution of temperature $T(x, t)$.

## Governing Equation
The 1D heat equation is given by:
$$\frac{\partial T}{\partial t} = \alpha \frac{\partial^2 T}{\partial x^2}$$
where:
- $T$: Temperature
- $t$: Time
- $x$: Spatial coordinate
- $\alpha$: Thermal diffusivity ($m^2/s$)

## Numerical Implementation
- **Spatial Discretization**: Second-order central difference for the Laplacian:
  $$\frac{\partial^2 T}{\partial x^2} \approx \frac{T_{i+1} - 2T_i + T_{i-1}}{\Delta x^2}$$
- **Time Integration**: Forward Euler (Explicit):
  $$T_i^{n+1} = T_i^n + \Delta t \cdot \alpha \left( \frac{T_{i+1}^n - 2T_i^n + T_{i-1}^n}{\Delta x^2} \right)$$
- **Stability**: Requires $\Delta t \le \frac{\Delta x^2}{2\alpha}$ (CFL condition).
