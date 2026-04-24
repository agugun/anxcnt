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
# # 1D Heat Equation (Implicit) - The Analytical Bridge
#
# This notebook demonstrates the complete workflow from mathematical conception to high-performance C++ execution for the **1D Heat Conduction** problem.
#
# ---

# %% [markdown]
# ## 1. Conceptual Framework
#
# The 1D Heat Equation describes how temperature $T$ evolves in space $x$ and time $t$ through a material with thermal diffusivity $\alpha$.
#
# ### Governing PDE
# $$\frac{\partial T}{\partial t} = \alpha \frac{\partial^2 T}{\partial x^2}$$
#
# ### Boundary Conditions (Dirichlet)
# - $T(0, t) = T_{left}$
# - $T(L, t) = T_{right}$

# %% [markdown]
# ## 2. Analytical Trace (Sympy)
#
# We use **Separation of Variables** to find the exact solution.

# %%
import sympy as sp
from IPython.display import display, Math

sp.init_printing()
x, t, alpha, L = sp.symbols('x t alpha L', real=True, positive=True)
T = sp.Function('T')(x, t)
X = sp.Function('X')(x)
Phi = sp.Function('Phi')(t)
lam = sp.symbols('lambda', real=True, positive=True)

spatial_ode = sp.Eq(X.diff(x, x) + lam**2 * X, 0)
temporal_ode = sp.Eq(Phi.diff(t) + alpha * lam**2 * Phi, 0)

display(Math(f"\\text{{Spatial ODE: }} {sp.latex(spatial_ode)}"))
display(Math(f"\\text{{Temporal ODE: }} {sp.latex(temporal_ode)}"))

sol_x = sp.dsolve(spatial_ode, X)
display(Math(f"X(x) = {sp.latex(sol_x.rhs)}"))

# %% [markdown]
# ## 3. Numerical Prototyping (Numpy/JAX)
#
# Before moving to C++, we prototype the **Implicit Euler** scheme in Python.

# %%
import numpy as np
import jax.numpy as jnp
from jax import jit
from scipy.linalg import solve_banded

# Parameters
nx = 100
L_val = 1.0
alpha_val = 0.05
dt = 0.01
dx = L_val / (nx - 1)
r = alpha_val * dt / dx**2

x_vals = np.linspace(0, L_val, nx)
u0 = np.exp(-100 * (x_vals - 0.5)**2)
u0[0] = u0[-1] = 0.0

def solve_implicit_numpy(u, nt):
    u_curr = u.copy()
    # Tridiagonal Matrix: (1+2r) on diag, -r on sub/super diag
    diag = (1 + 2*r) * np.ones(nx-2)
    off_diag = -r * np.ones(nx-3)
    ab = np.zeros((3, nx-2))
    ab[0, 1:] = off_diag
    ab[1, :] = diag
    ab[2, :-1] = off_diag
    
    for _ in range(nt):
        rhs = u_curr[1:-1]
        u_curr[1:-1] = solve_banded((1, 1), ab, rhs)
    return u_curr


# %% [markdown]
# ## 4. Performance & Verification Benchmark
#
# We compare our JAX/Numpy prototypes against calculations from the high-performance **axcnt_cpp** engine.

# %%
import sys, os, time
sys.path.append(os.path.abspath('../../../bindings/python'))
import axcnt_cpp
import matplotlib.pyplot as plt

# 1. Run C++ Backend
T_cond = [alpha_val] * (nx - 1)
storage = [1.0] * nx
cpp_sim = axcnt_cpp.Heat1D(T_cond, storage, 0.0, 0.0)
cpp_sim.set_initial_condition(list(u0))

start = time.time()
nt = 100
for _ in range(nt):
    cpp_sim.step(dt)
u_cpp = cpp_sim.get_values()
print(f"C++ Execution Time: {time.time() - start:.4f}s")

# 2. Run Numpy Reference
start = time.time()
u_np = solve_implicit_numpy(u0, nt)
print(f"Numpy Execution Time: {time.time() - start:.4f}s")

# 3. Verification Plot
plt.figure(figsize=(10, 6))
plt.plot(x_vals, u0, 'k:', label='Initial condition')
plt.plot(x_vals, u_np, 'b-', label='Numerical (Numpy)')
plt.plot(x_vals, u_cpp, 'ro', markersize=4, alpha=0.5, label='C++ Backend Parity')
plt.title("Analytical Bridge Parity: 1D Heat Equation")
plt.legend()
plt.grid(True)
plt.show()
