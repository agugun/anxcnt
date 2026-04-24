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
# # 2D Fluid Dynamics - Finite Element Analytical Bridge
#
# This notebook demonstrates the mathematical modeling and numerical simulation of **2D Incompressible Flow** using the Finite Element Method (FEM).
#
# ---

# %% [markdown]
# ## 1. Conceptual Framework
#
# Steady-state incompressible flow at low Reynolds numbers is governed by the **Stokes Equations**.
#
# ### Governing PDEs
# - **Momentum**: $-\mu \nabla^2 \mathbf{u} + \nabla p = \mathbf{f}$
# - **Continuity**: $\nabla \cdot \mathbf{u} = 0$
#
# Where $\mathbf{u}$ is the velocity vector and $p$ is the pressure field.

# %% [markdown]
# ## 2. Analytical Trace: Variational Form (Sympy)
#
# We derive the **Weak Form** used for Finite Element implementation.

# %%
import sympy as sp
from IPython.display import display, Math

mu = sp.symbols('mu', real=True, positive=True)
x, y = sp.symbols('x y')
u = sp.Function('u')(x, y)
v = sp.Function('v')(x, y) # Test function

# Laplacian (Diffusion) Term Weak Form
weak_diffusion = mu * (u.diff(x)*v.diff(x) + u.diff(y)*v.diff(y))
display(Math(f"\\int_{{\\Omega}} {sp.latex(weak_diffusion)} d\\Omega"))

# %% [markdown]
# ## 3. Numerical Prototyping (JAX)
#
# We prototype a simple stream-function solver to visualize flow patterns before full FEM execution.

# %%
import numpy as np
import jax.numpy as jnp
from jax import jit

nx, ny = 40, 40
dx, dy = 1.0/(nx-1), 1.0/(ny-1)

@jit
def sor_step(psi):
    # Solve Poisson eq for streamfunction: Del^2 psi = -omega
    psi_new = 0.25 * (jnp.roll(psi, 1, 0) + jnp.roll(psi, -1, 0) + 
                      jnp.roll(psi, 1, 1) + jnp.roll(psi, -1, 1))
    return psi_new


# %% [markdown]
# ## 4. Performance & Verification Benchmark
#
# Comparison with the **high-performance C++ FEM backend** (Stabilized P1-P1 ELEMENTS).

# %%
import sys, os, time
sys.path.append(os.path.abspath('../../bindings/python'))
import axcnt_cpp
import matplotlib.pyplot as plt

# 1. Initialize C++ Stokes Engine
nx_val, ny_val = 30, 30
cpp_sim = axcnt_cpp.Stokes2D(nx_val, ny_val, 1.0, 1.0, 0.1)

# 2. Set BCs (Lid Driven Cavity)
for i in range(nx_val):
    # Top boundary: u=1, v=0
    top_node = (ny_val-1)*nx_val + i
    cpp_sim.set_boundary_condition(top_node, 1.0, 0.0)
    # Others are default 0

start = time.time()
cpp_sim.solve()
print(f"C++ FEM Solve Time: {time.time() - start:.4f}s")

u_cpp = np.array(cpp_sim.get_u()).reshape(ny_val, nx_val)
v_cpp = np.array(cpp_sim.get_v()).reshape(ny_val, nx_val)

plt.figure(figsize=(8, 8))
X, Y = np.meshgrid(np.linspace(0, 1, nx_val), np.linspace(0, 1, ny_val))
plt.streamplot(X, Y, u_cpp, v_cpp, color='blue', density=1.5)
plt.title("Analytical Bridge: C++ Stokes FEM Streamlines (Lid-Driven)")
plt.show()
