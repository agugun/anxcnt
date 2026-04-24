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
# # 1D Wave Equation - Propagation Analytical Bridge
#
# This notebook demonstrates the modeling of mechanical wave propagation, transitioning from D'Alembert's analytical principles to C++ numerical execution.
#
# ---

# %% [markdown]
# ## 1. Conceptual Framework
#
# The 1D Wave Equation describes the displacement $u$ of a medium over time and space.
#
# ### Governing PDE
# $$\frac{\partial^2 u}{\partial t^2} = c^2 \frac{\partial^2 u}{\partial x^2}$$
#
# Where $c$ is the wave speed. This is a hyperbolic PDE, allowing for wave-front preservation.

# %% [markdown]
# ## 2. Analytical Trace: D'Alembert Solution (Sympy)
#
# We verify that any function of type $f(x - ct)$ or $g(x + ct)$ satisfies the wave equation.

# %%
import sympy as sp
from IPython.display import display, Math

x, t, c = sp.symbols('x t c', real=True, positive=True)
f = sp.Function('f')
u_ansatz = f(x - c*t)

lhs = u_ansatz.diff(t, t)
rhs = c**2 * u_ansatz.diff(x, x)

display(Math(f"\\frac{{\\partial^2 u}}{{\\partial t^2}} = {sp.latex(lhs)}"))
display(Math(f"c^2 \\frac{{\\partial^2 u}}{{\\partial x^2}} = {sp.latex(rhs)}"))
display(Math(f"\\text{{Parity: }} {sp.simplify(lhs - rhs) == 0}"))

# %% [markdown]
# ## 3. Numerical Prototyping (JAX)
#
# We prototype the second-order central difference scheme.

# %%
import jax.numpy as jnp
from jax import jit
import numpy as np

nx = 100
L = 1.0
dx = L / (nx - 1)
c_val = 1.0
dt = 0.5 * dx / c_val # CFL < 1
r2 = (c_val * dt / dx)**2

@jit
def wave_step(u, u_prev):
    return 2*u - u_prev + r2 * (jnp.roll(u, -1) - 2*u + jnp.roll(u, 1))


# %% [markdown]
# ## 4. Performance & Verification Benchmark
#
# Comparison with the **high-performance C++ backend** using the `Wave1D` bridge.

# %%
import sys, os, time
sys.path.append(os.path.abspath('../../bindings/python'))
import axcnt_cpp
import matplotlib.pyplot as plt

x_vals = np.linspace(0, L, nx)
u0 = np.exp(-100 * (x_vals - 0.3)**2)

# 1. Initialize C++ Engine
T_cond = [c_val**2/dx**2] * (nx - 1) 
storage = [1.0] * nx
cpp_sim = axcnt_cpp.Wave1D(T_cond, storage)
cpp_sim.set_initial_condition(list(u0), [0.0]*nx)

nt = 40
start = time.time()
for _ in range(nt):
    cpp_sim.step(dt)
u_cpp = cpp_sim.get_values()[:nx]
print(f"C++ Wave Execution Time: {time.time() - start:.4f}s")

plt.figure(figsize=(10, 6))
plt.plot(x_vals, u0, 'k:', label='Initial pulse')
plt.plot(x_vals, u_cpp, 'r-', label='C++ Backend (Propagated)')
plt.title("Wave Propagation: C++ Backend Parity")
plt.legend()
plt.grid(True)
plt.show()
