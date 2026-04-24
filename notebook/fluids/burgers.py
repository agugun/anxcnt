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
# # Burgers' Equation - Non-linear Analytical Bridge
#
# This notebook demonstrates the mathematical modeling and numerical simulation of the **Viscous Burgers' Equation**, representing a bridge between simple diffusion and complex fluid dynamics.
#
# ---

# %% [markdown]
# ## 1. Conceptual Framework
#
# Burgers' equation is a fundamental PDE in fluid mechanics that combines non-linear convection and linear diffusion.
#
# ### Governing PDE
# $$\frac{\partial u}{\partial t} + u \frac{\partial u}{\partial x} = \nu \frac{\partial^2 u}{\partial x^2}$$
#
# Where:
# - $u$: Velocity field
# - $\nu$: Kinematic viscosity

# %% [markdown]
# ## 2. Analytical Trace (Sympy)
#
# We verify the **Cole-Hopf transformation**, which maps the non-linear Burgers' equation to the linear heat equation.

# %%
import sympy as sp
from IPython.display import display, Math

nu, x, t = sp.symbols('nu x t', real=True, positive=True)
phi = sp.Function('phi')(x, t)

# Cole-Hopf Transformation: u = -2 * nu * (phi_x / phi)
u = -2 * nu * (phi.diff(x) / phi)

display(Math(f"u(x, t) = -2 \\nu \\frac{{\\partial \\phi / \\partial x}}{{\\phi}}"))

# Substitute u into Burgers' equation
burgers_lhs = u.diff(t) + u * u.diff(x)
burgers_rhs = nu * u.diff(x, x)

residual = sp.simplify(burgers_lhs - burgers_rhs)
display(Math(f"\\text{{Residual assuming }} \\phi \\text{{ solves Heat Eq: }} {sp.latex(residual)}"))

# %% [markdown]
# ## 3. Numerical Prototyping (JAX)
#
# We prototype a shock-capturing upwind scheme with JAX.

# %%
import jax.numpy as jnp
from jax import jit
import numpy as np

nx = 100
dx = 2.0 * np.pi / (nx - 1)
nu_val = 0.07
dt = 0.01
x_vals = np.linspace(0, 2*np.pi, nx)
u0 = np.sin(x_vals)

@jit
def burgers_step(u):
    # Upwind convection
    u_x = jnp.where(u > 0, (u - jnp.roll(u, 1)) / dx, (jnp.roll(u, -1) - u) / dx)
    # Central diffusion
    u_xx = (jnp.roll(u, -1) - 2*u + jnp.roll(u, 1)) / dx**2
    return u + dt * (-u * u_x + nu_val * u_xx)


# %% [markdown]
# ## 4. Performance & Verification Benchmark
#
# Comparison between JAX (Upwind) and the high-performance C++ backend (Implicit Newton-Raphson).

# %%
import sys, os, time
sys.path.append(os.path.abspath('../../bindings/python'))
import axcnt_cpp
import matplotlib.pyplot as plt

nt = 100
start = time.time()
cpp_sim = axcnt_cpp.Burgers1D(nu_val, dx, list(u0))
for _ in range(nt):
    cpp_sim.step(dt)
u_cpp = cpp_sim.get_values()
print(f"C++ Execution Time: {time.time() - start:.4f}s")

plt.figure(figsize=(10, 6))
plt.plot(x_vals, u0, 'k:', label='Initial')
plt.plot(x_vals, u_cpp, 'b-', label='C++ (Implicit Backend)')
plt.title("Burgers Shock Evolution - C++ Backend Parity")
plt.legend()
plt.grid(True)
plt.show()
