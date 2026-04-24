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
# # 2D Heat Equation (Implicit) - Multi-Dimensional Bridge
#
# This notebook demonstrates the mathematical modeling and numerical simulation of **2D Heat Conduction** in a rectangular domain.
#
# ---

# %% [markdown]
# ## 1. Conceptual Framework
#
# The 2D Heat Equation is an extension of the 1D model into two spatial dimensions:
# $$\frac{\partial T}{\partial t} = \alpha \left( \frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2} \right)$$
#
# ### Boundary Conditions (Dirichlet)
# We assume fixed temperatures at all four boundaries: $T_{top}, T_{bottom}, T_{left}, T_{right}$.

# %% [markdown]
# ## 2. Analytical Trace (Sympy)
#
# We derive the 2D fundamental harmonic solution.

# %%
import sympy as sp
from IPython.display import display, Math

sp.init_printing()
x, y, t, alpha, Lx, Ly = sp.symbols('x y t alpha Lx Ly', real=True, positive=True)
n, m = sp.symbols('n m', integer=True, positive=True)

# Fundamental harmonic
T_nm = sp.sin(n * sp.pi * x / Lx) * sp.sin(m * sp.pi * y / Ly) * sp.exp(-alpha * ((n*sp.pi/Lx)**2 + (m*sp.pi/Ly)**2) * t)

display(Math(f"T_{{n,m}}(x, y, t) = A_{{n,m}} {sp.latex(T_nm)}"))

# %% [markdown]
# ## 3. Numerical Prototyping (JAX)
#
# 2D problems significantly increase the number of grid points ($N \times M$). We use JAX for vectorized prototyping.

# %%
import jax.numpy as jnp
from jax import jit
import numpy as np

nx, ny = 50, 50
Lx, Ly = 1.0, 1.0
dx, dy = Lx/(nx-1), Ly/(ny-1)
alpha_val = 0.05
dt = 0.001

@jit
def heat_2d_step(T):
    d2T_dx2 = (jnp.roll(T, -1, axis=1) - 2*T + jnp.roll(T, 1, axis=1)) / dx**2
    d2T_dy2 = (jnp.roll(T, -1, axis=0) - 2*T + jnp.roll(T, 1, axis=0)) / dy**2
    return T + alpha_val * dt * (d2T_dx2 + d2T_dy2)


# %% [markdown]
# ## 4. Performance & Verification Benchmark
#
# We compare the Python prototype results with the C++ backend's 2D engine.

# %%
import sys, os, time
sys.path.append(os.path.abspath('../../../bindings/python'))
import axcnt_cpp
import matplotlib.pyplot as plt

# 1. Initialize C++ Backend
cpp_sim = axcnt_cpp.Heat2D(nx, ny, Lx, Ly, alpha_val)
u0 = np.zeros((ny, nx))
u0[ny//4:3*ny//4, nx//4:3*nx//4] = 1.0
cpp_sim.set_initial_condition(list(u0.flatten()))

# 2. Run C++ Simulation
start = time.time()
nt = 100
for _ in range(nt):
    cpp_sim.step(dt)
u_cpp = np.array(cpp_sim.get_values()).reshape(ny, nx)
print(f"C++ Execution Time: {time.time() - start:.4f}s")

# 3. Plot Result
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.imshow(u0, cmap='hot', origin='lower')
plt.title("Initial Condition")
plt.subplot(1, 2, 2)
plt.imshow(u_cpp, cmap='hot', origin='lower')
plt.title(f"C++ Result after {nt} steps")
plt.show()
