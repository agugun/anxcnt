# 2D Dual-Phase (Oil-Water) Reservoir Simulation

## Scenario
This module simulates a 2D waterflooding scenario. A water injector is placed at one corner of the reservoir, and an oil producer is at the opposite corner. The simulation tracks the pressure distribution and the movement of the water saturation front as it sweeps oil toward the producer.

## Mathematics
The simulation uses a **Fully Implicit (FI)** formulation for coupled Oil-Water flow:

1.  **Non-Linear Residuals**:
    - **Water**: $R_{w,i} = \frac{\phi V}{\Delta t} (S_{w,i}^{n+1} - S_{w,i}^n) - \sum T_{w,ij} (P_j - P_i) - Q_{w,i} = 0$
    - **Oil**: $R_{o,i} = \frac{\phi V}{\Delta t} (1 - S_{w,i}^{n+1} - (1 - S_{w,i}^n)) - \sum T_{o,ij} (P_j - P_i) - Q_{o,i} = 0$

2.  **Newton-Raphson Solver**:
    Iteratively solves for $X = [P, S_w]$ using the $2 \times 2$ block Jacobian:
    $$X_{k+1} = X_k - J^{-1} R(X_k)$$

### Relative Permeability
Uses a **Corey-type Model**:
- $k_{rw} = S_{we}^{2}$
- $k_{ro} = (1-S_{we})^{2}$

## Logic
- **Fully Implicit**: Both pressure and saturation are solved simultaneously. This provides unconditional stability (up to Newton convergence) and handles the coupling between $P$ and $S_w$ exactly.
- **Newton Solver**: Implemented with **Step Length Damping** to handle the non-linear displacement front, especially during start-up.
- **Upwinding**: Local transmissibilities for both phases use the saturation of the higher-pressure (upwind) block.
- **Visualization**: Generates multi-field `.vti` files containing both Pressure and WaterSaturation.
- **Stability**: Much more robust than IMPES; handles time steps of **12-24 hours** without numerical oscillations.
