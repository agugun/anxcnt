# 2D Dual-Phase Oil and Gas Simulation

This module simulates the displacement of viscous oil by a highlight compressible gas phase in a 2D porous medium. It captures the complex physics of primary depletion and mixed-phase flow.

## Scenario Description

This module simulates **Primary Depletion** in a 2D reservoir initially containing both oil and gas.

- **Grid**: 20x20 ($2000 \times 2000$ ft)
- **Thickness**: 50 ft
- **Initial Conditions**: 
    - Uniform pressure: 5000 psi
    - Uniform gas saturation: $S_{gi} = 0.15$ (with $S_{wc} = 0.2$ fixed)
- **Wells**:
    - **Single Producer**: Located at (19,19), producing 200 STB/D of reservoir fluid.

## Mathematical Formulation

The simulation uses a **Fully Implicit (Black Oil Lite)** approach to solve the coupled non-linear mass balance equations:

### 1. Oil Mass Balance
$$ \nabla \cdot \left( \frac{k k_{rog}}{\mu_o B_o} \nabla P_o \right) = \frac{\partial}{\partial t} \left( \frac{\phi S_o}{B_o} \right) + q_o $$

### 2. Gas Mass Balance
$$ \nabla \cdot \left( \frac{k k_{rg}}{\mu_g B_g} \nabla P_o \right) = \frac{\partial}{\partial t} \left( \frac{\phi S_g}{B_g} \right) + q_g $$

### Key Property Models
- **Gas Compressibility ($B_g$)**: Modeled via the Ideal Gas Law: $B_g \propto 1/P$.
- **Relative Permeability**: Gas-Oil Corey curves with $S_{wc}=0.2$ and $S_{org}=0.1$.
- **Upwinding**: Integrated upstream weighting for flux terms to ensure stability.

## Performance Notes
Due to the high compressibility of gas, the system uses the `NewtonSolver` with a direct `LUSolver` to ensure robust convergence during pressure depletion.
