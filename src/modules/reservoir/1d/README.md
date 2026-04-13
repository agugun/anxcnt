# 1D Single-Phase Reservoir Simulation

## Scenario
This module simulates a 1D linear reservoir depletion. A single production well is located at the center of a linear strip (e.g., a core plug or a narrow reservoir section). The reservoir is initially at a higher pressure and is produced at a constant surface rate.

## Mathematics
The simulation solves the **Pressure Diffusivity Equation** in 1D:
$$\frac{\partial p}{\partial t} = \eta \frac{\partial^2 p}{\partial x^2} - \frac{Q B}{\phi c_t V}$$

Where:
- $\eta = 0.0002637 \frac{k}{\phi \mu c_t}$ (Diffusivity in $ft^2/hr$)
- $k$ is permeability [mD]
- $\phi$ is porosity [fraction]
- $\mu$ is viscosity [cP]
- $c_t$ is total compressibility [$psi^{-1}$]

## Logic
- **Integration**: Uses the **Implicit Euler** method for unconditional stability.
- **Solver**: Implements the **Thomas Algorithm** (Linear Tridiagonal Solver), which is $O(N)$ and optimal for 1D systems.
- **Boundaries**: **No-flow (Neumann)** conditions are applied at both ends ($\frac{\partial p}{\partial x} = 0$), simulating a closed reservoir.
- **Units**: Standard Oil Field Units.
