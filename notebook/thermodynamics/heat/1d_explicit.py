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
# # 1D Heat Equation (Explicit) - Analytical to Numerical Bridge
#
# This notebook explores the **Explicit Forward Euler** scheme for the heat equation, tracing the path from symbolic derivation to high-performance execution.
#
# ---

# %% [markdown]
# ## 1. Conceptual Framework
#
# The 1D Heat Equation:
# $$\frac{\partial T}{\partial t} = \alpha \frac{\partial^2 T}{\partial x^2}$$
#
# ### Numerical Discretization (Explicit)
# Using Forward Euler in time and Central Difference in space:
# $$T_i^{n+1} = T_i^n + \Delta t \cdot \alpha \left( \frac{T_{i+1}^n - 2T_i^n + T_{i-1}^n}{\Delta x^2} \right)$$
#
# **Stability Condition**: $Fo = \frac{\alpha \Delta t}{\Delta x^2} \le \frac{1}{2}$.

# %% [markdown]
# ## 2. Analytical Trace (Sympy)
#
# We verify the fundamental solution $T(x, t) = A \sin(kx) e^{-\alpha k^2 t}$.

# %%
import sympy as sp
from IPython.display import display, Math

x, t, alpha, k, A = sp.symbols('x t alpha k A', real=True)
T_expr = A * sp.sin(k*x) * sp.exp(-alpha * k**2 * t)

# Verify PDE satisfaction
lhs = T_expr.diff(t)
rhs = alpha * T_expr.diff(x, x)
display(Math(f"{sp.latex(lhs)} = {sp.latex(rhs)}"))
display(Math(f"\\text{{Residual: }} {sp.simplify(lhs - rhs)}"))

# %% [markdown]
# ## 3. Numerical Prototyping (JAX)
#
# We use JAX for high-performance Python prototyping before moving to C++.

# %%
import numpy as np
import jax.numpy as jnp
from jax import jit

nx = 100
L = 1.0
dx = L / (nx - 1)
alpha_val = 0.05
dt = 0.4 * dx**2 / alpha_val # Safe for explicit
r = alpha_val * dt / dx**2

x_vals = np.linspace(0, L, nx)
u0 = np.exp(-100 * (x_vals - 0.5)**2)

@jit
def explicit_step(u):
    # Standard central difference stencil
    return u + r * (jnp.roll(u, -1) - 2*u + jnp.roll(u, 1))

def jax_solver(u, steps):
    u_curr = jnp.array(u)
    for _ in range(steps):
        u_curr = explicit_step(u_curr)
        # Dirichlet Boundary
        u_curr = u_curr.at[0].set(0.0)
        u_curr = u_curr.at[-1].set(0.0)
    return u_curr


# %% [markdown]
# ## 4. Performance & Verification Benchmark
#
# Comparison between JAX (Explicit) and the high-performance C++ backend (Implicit).

# %%
import sys, os, time
sys.path.append(os.path.abspath('../../../bindings/python'))
import axcnt_cpp
import matplotlib.pyplot as plt

nt = 100
start = time.time()
u_jax = jax_solver(u0, nt).block_until_ready()
print(f"JAX Execution Time: {time.time() - start:.4f}s")

T_cond = [alpha_val] * (nx - 1)
cpp_sim = axcnt_cpp.Heat1D(T_cond, [1.0]*nx, 0.0, 0.0)
cpp_sim.set_initial_condition(list(u0))

start = time.time()
for _ in range(nt):
    cpp_sim.step(dt)
u_cpp = cpp_sim.get_values()
print(f"C++ Execution Time: {time.time() - start:.4f}s")

plt.figure(figsize=(10, 6))
plt.plot(x_vals, u0, 'k:', label='Initial')
plt.plot(x_vals, u_jax, 'b-', label='JAX (Explicit)')
plt.plot(x_vals, u_cpp, 'ro', markersize=4, alpha=0.5, label='C++ (Implicit Backend)')
plt.title("Explicit vs Implicit Parity")
plt.legend()
plt.show()
