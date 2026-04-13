# Pressure Diffusivity 1D

## Physics Background
This module simulates the 1D pressure transient in a porous medium (reservoir). It models how fluid pressure evolves in a reservoir due to production or injection at boundaries.

## Governing Equation
The 1D pressure diffusivity equation for an undersaturated liquid is:
$$\frac{\partial P}{\partial t} = \eta \frac{\partial^2 P}{\partial x^2}$$
where $\eta$ is the diffusivity constant:
$$\eta = 0.0002637 \frac{k}{\phi \mu c_t} \text{ (Field Units)}$$
- $k$: Permeability (md)
- $\phi$: Porosity (fraction)
- $\mu$: Viscosity (cp)
- $c_t$: Total compressibility (1/psi)

## Numerical Implementation
- **Scheme**: Implicit Euler (Forward in Time, Central in Space).
- **Solver**: **Linear Tridiagonal Solver** (Thomas Algorithm).
- **Efficiency**: The implicit formulation allows for large reservoir-scale time steps (days/weeks) while maintaining numerical stability.
