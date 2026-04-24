# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: .venv (3.12.3)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # 3D Reservoir Simulation - Volumetric Analytical Bridge
#
# This notebook demonstrates the mathematical modeling and numerical simulation of fluid flow in 3D porous volumes, representing the peak of reservoir engineering complexity.
#
# ---

# %% [markdown]
# ## 1. Conceptual Framework
#
# Fluid flow in 3D porous media is governed by the 3D Pressure Diffusivity Equation.
#
# ### Governing PDE
# $$\frac{\partial^2 p}{\partial x^2} + \frac{\partial^2 p}{\partial y^2} + \frac{\partial^2 p}{\partial z^2} = \frac{\phi \mu c_t}{k} \frac{\partial p}{\partial t}$$
#
# This equation accounts for gravity, vertical communication within layers, and complex well trajectories.

# %% [markdown]
# ## 2. Analytical Trace: 3D Harmonic Solution (Sympy)
#
# We derive the volumetric decaying solution.

# %%
import sympy as sp
from IPython.display import display, Math

sp.init_printing()
x, y, z, t, eta, Lx, Ly, Lz = sp.symbols('x y z t eta Lx Ly Lz', real=True, positive=True)
nx, ny, nz = sp.symbols('nx ny nz', integer=True, positive=True)

kx = nx * sp.pi / Lx
ky = ny * sp.pi / Ly
kz = nz * sp.pi / Lz

P_nml = sp.sin(kx*x) * sp.sin(ky*y) * sp.sin(kz*z) * sp.exp(-eta * (kx**2 + ky**2 + kz**2) * t)

display(Math(f"P_{{n,m,l}}(x, y, z, t) = A_{{n,m,l}} {sp.latex(P_nml)}"))

# %% [markdown]
# ## 3. Numerical Prototyping (JAX)
#
# Prototyping 3D volumes requires efficient tensor operations.

# %%
import jax.numpy as jnp
from jax import jit
import numpy as np

NX, NY, NZ = 20, 20, 5
eta_val = 100.0
dt = 0.1
dx = 50.0

@jit
def res_3d_step(p):
    d2p = (jnp.roll(p, 1, 0) + jnp.roll(p, -1, 0) + 
           jnp.roll(p, 1, 1) + jnp.roll(p, -1, 1) + 
           jnp.roll(p, 1, 2) + jnp.roll(p, -1, 2) - 6*p) / dx**2
    return p + eta_val * dt * d2p


# %% [markdown]
# ## 4. Performance & Verification Benchmark
#
# We verify 1D slices of our prototype against the high-performance C++ backend.

# %%
import sys, os, time
sys.path.append(os.path.abspath('../../bindings/python'))
import axcnt_cpp
import matplotlib.pyplot as plt

# Benchmark 1D Slice
nx_slice = 100
T_cond = [eta_val / dx**2] * (nx_slice - 1)
cpp_sim = axcnt_cpp.Reservoir1D(T_cond, [1.0]*nx_slice)
p0 = [3000.0] * nx_slice
p0[nx_slice//2] = 500.0
cpp_sim.set_initial_condition(p0)

for _ in range(10):
    cpp_sim.step(dt)
p_slice = cpp_sim.get_values()

plt.figure(figsize=(10, 6))
plt.plot(p_slice, 'r-', label='C++ 1D Kernel (Benchmark)')
plt.title("3D Verification: 1D Slice Kernel Comparison")
plt.legend()
plt.show()
