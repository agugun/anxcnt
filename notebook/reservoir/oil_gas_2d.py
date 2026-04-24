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
# # Oil-Gas Reservoir Simulation - Multi-phase Analytical Bridge
#
# This notebook demonstrates the mathematical modeling and numerical simulation of **Two-phase Flow** (Oil-Gas) in porous media.
#
# ---

# %% [markdown]
# ## 1. Conceptual Framework
#
# In multi-phase flow, we must track the saturation $S_j$ and pressure $p_j$ for each phase $j$.
#
# ### Accumulation Law
# For the oil phase:
# $$\frac{\partial}{\partial t} \left( \frac{\phi S_o}{B_o} \right) = \nabla \cdot \left( \frac{k k_{ro}}{\mu_o B_o} \nabla p_o \right)$$
#
# For the gas phase:
# $$\frac{\partial}{\partial t} \left( \frac{\phi S_g}{B_g} + \frac{\phi R_s S_o}{B_o} \right) = \nabla \cdot \left( \dots \right)$$
#
# This system is highly non-linear due to relative permeability $k_{rj}(S_j)$.

# %% [markdown]
# ## 2. Analytical Trace: Expansion Terms (Sympy)
#
# We use Sympy to expand the temporal accumulation derivatives into pressure and saturation components.

# %%
import sympy as sp
from IPython.display import display, Math

t, phi, So, Bo, Sg, Bg, Rs = sp.symbols('t phi S_o B_o S_g B_g R_s', real=True)

# Accumulation component for Oil
acc_o = phi * So / Bo
display(Math(f"\\text{{Oil Accumulation: }} {sp.latex(acc_o)}"))

# Expansion using chain rule
display(Math(f"\\frac{{\\partial}}{{\\partial t}}() = {sp.latex(sp.diff(acc_o, t))}"))

# %% [markdown]
# ## 3. Numerical Prototyping (JAX)
#
# Prototyping phase transition logic with JAX.

# %%
import jax.numpy as jnp
from jax import jit

@jit
def bubble_point_check(p, rs_vals, pb_const):
    # Determine if free gas exists
    return jnp.where(p < pb_const, "Undefined Saturation Change", "Undersaturated")


# %% [markdown]
# ## 4. Performance & Verification Benchmark
#
# Comparison with the **high-performance C++ backend (2D Pressure Engine)** for the pressure propagation component.

# %%
import sys, os, time
sys.path.append(os.path.abspath('../../bindings/python'))
import axcnt_cpp
import numpy as np
import matplotlib.pyplot as plt

nx, ny = 30, 30
cpp_sim = axcnt_cpp.Reservoir2D(nx, ny, 3000.0, 3000.0, 100.0)
cpp_sim.set_initial_condition([3000.0] * (nx*ny))

start = time.time()
cpp_sim.step(1.0)
p_res = cpp_sim.get_values()
print(f"C++ Oil-Gas Pressure Kernel Solve Time: {time.time() - start:.4f}s")

plt.figure(figsize=(10, 6))
plt.contourf(np.array(p_res).reshape(ny, nx))
plt.title("Analytical Bridge: C++ Multi-phase Pressure Kernel")
plt.show()
