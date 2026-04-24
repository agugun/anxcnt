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
# # 1D Reservoir Simulation - Single Phase Analytical Bridge
#
# This notebook demonstrates the mathematical modeling and numerical simulation of fluid flow in a 1D reservoir, bridging the gap from simple diffusion to engineered solutions.
#
# ---

# %% [markdown]
# ## 1. Conceptual Framework
#
# Reservoir flow in 1D is governed by the conservation of mass and Darcy's Law.
#
# ### Governing PDE
# $$\frac{\partial}{\partial x} \left( \frac{k A}{\mu B} \frac{\partial p}{\partial x} \right) = A \phi c_t \frac{\partial p}{\partial t}$$
#
# In its simplest form, this recovers the diffusivity equation.

# %% [markdown]
# ## 2. Analytical Trace (Sympy)
#
# We use Sympy to solve the steady-state case (Laplace Equation) for a constant pressure boundary.

# %%
import sympy as sp
from IPython.display import display, Math

x, L, P1, P2 = sp.symbols('x L P1 P2', real=True)
p = sp.Function('p')(x)

steady_state_eq = sp.Eq(p.diff(x, x), 0)
sol = sp.dsolve(steady_state_eq, p, ics={p.subs(x, 0): P1, p.subs(x, L): P2})

display(Math(f"\\text{{Steady State Profile: }} {sp.latex(sol)}"))

# %% [markdown]
# ## 3. Numerical Prototyping (JAX)
#
# Prototyping a transient reservoir solve in Python.

# %%
import jax.numpy as jnp
from jax import jit
import numpy as np

nx = 50
dx = 100.0 # ft
x_vals = np.linspace(0, nx*dx, nx)
p0 = np.full(nx, 4000.0)

# %% [markdown]
# ## 4. Performance & Verification Benchmark
#
# Comparison with the **high-performance C++ backend** using the `Reservoir1D` bridge.

# %%
import sys, os, time
sys.path.append(os.path.abspath('../../bindings/python'))
import axcnt_cpp
import matplotlib.pyplot as plt

# 1. Initialize C++ Engine
T_cond = [5.0] * (nx - 1)  # trans
storage = [0.1] * nx # stor
cpp_sim = axcnt_cpp.Reservoir1D(T_cond, storage)
cpp_sim.set_initial_condition(list(p0))

nt = 50
dt = 1.0
start = time.time()
for _ in range(nt):
    cpp_sim.step(dt)
p_cpp = cpp_sim.get_values()
print(f"C++ Reservoir Execution Time: {time.time() - start:.4f}s")

plt.figure(figsize=(10, 6))
plt.plot(x_vals, p0, 'k:', label='Initial')
plt.plot(x_vals, p_cpp, 'r-', label='C++ Backend')
plt.title("Reservoir 1D: C++ Backend Parity")
plt.ylabel("Pressure (psi)")
plt.legend()
plt.show()
