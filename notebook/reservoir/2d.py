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
# # 2D Reservoir Simulation - Spatial Analytical Bridge
#
# This notebook demonstrates the mathematical modeling and numerical simulation of fluid flow in a 2D reservoir domain.
#
# ---

# %% [markdown]
# ## 1. Conceptual Framework
#
# The 2D Reservoir Diffusivity equation describes pressure distribution $p(x, y, t)$ in a flat porous layer.
#
# ### Governing PDE
# $$\frac{\partial^2 p}{\partial x^2} + \frac{\partial^2 p}{\partial y^2} = \frac{\phi \mu c_t}{k} \frac{\partial p}{\partial t}$$
#
# This parabolic PDE captures the propagation of pressure waves from producers and injectors across the field.

# %% [markdown]
# ## 2. Analytical Trace: Separation of Variables (Sympy)
#
# We derive the eigenfunctions of a rectangular reservoir with constant pressure boundaries.

# %%
import sympy as sp
from IPython.display import display, Math

sp.init_printing()
x, y, t, eta, Lx, Ly = sp.symbols('x y t eta Lx Ly', real=True, positive=True)
n, m = sp.symbols('n m', integer=True, positive=True)

# Fundamental harmonic for 2D Pressure Diffusion
p_nm = sp.sin(n * sp.pi * x / Lx) * sp.sin(m * sp.pi * y / Ly) * sp.exp(-eta * ((n*sp.pi/Lx)**2 + (m*sp.pi/Ly)**2) * t)

display(Math(f"p_{{n,m}}(x, y, t) = A_{{n,m}} {sp.latex(p_nm)}"))

# %% [markdown]
# ## 3. Numerical Prototyping (JAX)
#
# Vectorized pressure prototyping with JAX.

# %%
import jax.numpy as jnp
from jax import jit
import numpy as np

nx, ny = 60, 60
eta_val = 200.0
dt = 0.05
dx, dy = 100.0, 100.0

@jit
def res_2d_step(p):
    d2p_dx2 = (jnp.roll(p, -1, 1) - 2*p + jnp.roll(p, 1, 1)) / dx**2
    d2p_dy2 = (jnp.roll(p, -1, 0) - 2*p + jnp.roll(p, 1, 0)) / dy**2
    return p + eta_val * dt * (d2p_dx2 + d2p_dy2)


# %% [markdown]
# ## 4. Performance & Verification Benchmark
#
# Comparison with the **high-performance C++ backend** using the `Reservoir2D` bridge.

# %%
import sys, os, time
sys.path.append(os.path.abspath('../../bindings/python'))
import axcnt_cpp
import matplotlib.pyplot as plt

# 1. Initialize C++ Engine
cpp_sim = axcnt_cpp.Reservoir2D(nx, ny, nx*dx, ny*dy, eta_val)
p_init = np.full((ny, nx), 3000.0)
p_init[ny//2, nx//2] = 500.0 # Single producer
cpp_sim.set_initial_condition(list(p_init.flatten()))

nt = 20
start = time.time()
for _ in range(nt):
    cpp_sim.step(dt)
p_cpp = np.array(cpp_sim.get_values()).reshape(ny, nx)
print(f"C++ 2D Reservoir Execution Time: {time.time() - start:.4f}s")

plt.figure(figsize=(10, 8))
plt.contourf(p_cpp, levels=20, cmap='RdYlBu')
plt.colorbar(label='Pressure (psi)')
plt.title("Analytical Bridge: C++ 2D Reservoir Pressure Field")
plt.show()
