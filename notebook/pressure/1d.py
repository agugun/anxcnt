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
# # 1D Pressure Diffusivity - Reservoir Analytical Bridge
#
# This notebook demonstrates the mathematical modeling and numerical simulation of fluid pressure evolution in a porous medium.
#
# ---

# %% [markdown]
# ## 1. Conceptual Framework
#
# Pressure evolution in porous media combines Darcy's Law and mass conservation.
#
# ### Governing PDE
# $$\frac{\partial^2 p}{\partial x^2} = \frac{\phi \mu c_t}{k} \frac{\partial p}{\partial t}$$
#
# This is a diffusion-type PDE where the 'pressure diffusivity' is $\eta = \frac{k}{\phi \mu c_t}$.

# %% [markdown]
# ## 2. Analytical Trace: Fundamental Solution (Sympy)
#
# We derive the Gaussian pressure decay from an initial impulse.

# %%
import sympy as sp
from IPython.display import display, Math

x, t, eta = sp.symbols('x t eta', real=True, positive=True)
p_fundamental = (1 / sp.sqrt(4 * sp.pi * eta * t)) * sp.exp(-x**2 / (4 * eta * t))

display(Math(f"p(x, t) = {sp.latex(p_fundamental)}"))

# Verification
diff_eq = p_fundamental.diff(t) - eta * p_fundamental.diff(x, x)
display(Math(f"\\text{{PDE Residual: }} {sp.simplify(diff_eq)}"))

# %% [markdown]
# ## 3. Numerical Prototyping (JAX)
#
# Prototyping the implicit pressure solution in Python.

# %%
import jax.numpy as jnp
from jax import jit
import numpy as np

nx = 100
L = 2000.0
dx = L / (nx - 1)
eta_val = 500.0
dt = 1.0

x_vals = np.linspace(0, L, nx)
p0 = np.full(nx, 3000.0)
p0[nx//2] = 1000.0 # Well drawdown

# %% [markdown]
# ## 4. Performance & Verification Benchmark
#
# Comparison with the **high-performance C++ backend** using the `Pressure1D` bridge.

# %%
import sys, os, time
sys.path.append(os.path.abspath('../../bindings/python'))
import axcnt_cpp
import matplotlib.pyplot as plt

# 1. Initialize C++ Engine
T_cond = [eta_val / dx**2] * (nx - 1)
storage = [1.0] * nx
cpp_sim = axcnt_cpp.Pressure1D(T_cond, storage, 3000.0, 3000.0)
cpp_sim.set_initial_condition(list(p0))

nt = 20
start = time.time()
for _ in range(nt):
    cpp_sim.step(dt)
p_cpp = cpp_sim.get_values()
print(f"C++ Pressure Execution Time: {time.time() - start:.4f}s")

plt.figure(figsize=(10, 6))
plt.plot(x_vals, p0, 'k:', label='Initial Pressure')
plt.plot(x_vals, p_cpp, 'r-', label='C++ Backend (Diffused)')
plt.title("Pressure Diffusivity: C++ Backend Parity")
plt.ylabel("Pressure (psi)")
plt.legend()
plt.grid(True)
plt.show()
