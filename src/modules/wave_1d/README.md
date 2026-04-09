# Wave Equation 1D

## Physics Background
This module simulates the propagation of a pressure or displacement wave in a 1D medium. Unlike diffusion, waves are non-dissipative (ideally) and preserve energy and shape.

## Governing Equation
The 1D wave equation is a second-order PDE:
$$\frac{\partial^2 u}{\partial t^2} = c^2 \frac{\partial^2 u}{\partial x^2}$$
where:
- $u$: Displacement or pressure
- $c$: Wave speed

To fit the first-order system framework, it is rewritten as:
$$\frac{\partial u}{\partial t} = v$$
$$\frac{\partial v}{\partial t} = c^2 \frac{\partial^2 u}{\partial x^2}$$

## Numerical Implementation
- **Scheme**: Runge-Kutta 4 (RK4).
- **Why RK4?**: Wave equations are oscillatory. Explicit Euler methods are often unstable for waves, whereas RK4 provides the high-order accuracy and stability needed to preserve the wave shape over time.
- **Stability**: Requires $c \frac{\Delta t}{\Delta x} \le 1.0$ (CFL condition).
