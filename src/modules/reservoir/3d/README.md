# 3D Single-Phase Reservoir Simulation

## Scenario
Models a full volumetric reservoir block with a producer located in the center of the 3D volume. This allows for studying vertical pressure gradients and gravity-independent drainage bulbs.

## Mathematics
3D Pressure Diffusivity:
$$\frac{\partial p}{\partial t} = \eta \nabla^2 p - \text{Wells}$$

Discretized using a **7-point stencil** (Central + 6 neighbors).

## Logic
- **Efficient Computations**: Leverages a matrix-free **Conjugate Gradient** solver. This avoids the memory bottleneck ($O(N^2)$) of storing large 3D matrices.
- **Boundaries**: **No-flow** conditions on all 6 faces (Top, Bottom, and all 4 sides).
- **Visualization**: The `.vti` output is designed for 3D tools like ParaView. Users can use **Clip** or **Volume Rendering** filters to inspect the internal pressure distribution.
- **Parameters**: Permeability and porosity are volumetric averages.
                        
