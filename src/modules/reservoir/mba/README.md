# Material Balance (MBA)

## Physics Background
This module simulates the pressure depletion of a reservoir modeled as a single tank (0D model). It is a fundamental tool in reservoir engineering used to estimate STOIIP (Stock Tank Oil Initially In Place) and predict future performance.

## Governing Equation
For an undersaturated oil reservoir, the relationship between pressure depletion and production is governed by the total compressibility of the system:
$$\Delta V = V \cdot c_t \cdot \Delta P$$
Expressed as a time-dependent ODE:
$$\frac{\partial P}{\partial t} = -\frac{q(t)}{V \cdot c_t}$$
where:
- $P$: Reservoir Pressure (psi)
- $q(t)$: Production rate (STB/day)
- $V$: Reservoir volume (equivalent to $N \cdot B_o$)
- $c_t$: Total compressibility (1/psi)

## Numerical Implementation
- **Scheme**: Forward Euler (Explicit).
- **Complexity**: $O(1)$ per time step.
- **Usage**: Provides a rapid first-order estimate of reservoir pressure decline under constant or variable production rates.
