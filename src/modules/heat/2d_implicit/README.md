# Heat Equation 2D (Implicit)

## Physics Background
This module simulates 2D heat conduction on a rectangular grid. It demonstrates how to handle higher-dimensional PDEs and large linear systems using iterative solvers.

## Governing Equation
$$\frac{\partial T}{\partial t} = \alpha \left( \frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2} \right)$$

## Numerical Implementation
- **Spatial Discretization**: 5-point stencil for the 2D Laplacian:
  $$\Delta T \approx \frac{T_{i+1,j} - 2T_{i,j} + T_{i-1,j}}{\Delta x^2} + \frac{T_{i,j+1} - 2T_{i,j} + T_{i,j-1}}{\Delta y^2}$$
- **Time Integration**: Implicit Euler.
- **Solver**: **Conjugate Gradient (CG)**. 
- **Matrix-Free Approach**: For efficiency, the Jacobian is not explicitly stored as a matrix. Instead, the `apply_jacobian` method computes the 5-point stencil on the fly, reducing memory overhead from $O(N^4)$ to $O(N^2)$.
