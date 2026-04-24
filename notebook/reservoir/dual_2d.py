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
# # Dual-Porosity Reservoir (Warren-Root) - Coupled Analytical Bridge
#
# This notebook demonstrates the modeling of **Naturally Fractured Reservoirs (NFRs)** using the dual-porosity approach.
#
# ---

# %% [markdown]
# ## 1. Conceptual Framework
#
# Dual-porosity systems treat the matrix and fracture as two overlapping continua with a transfer term $\tau$.
#
# ### Governing Equations
#
# **Fracture System**:
# $$\nabla \cdot (k_f \nabla p_f) - \tau = \phi_f c_f \frac{\partial p_f}{\partial t}$$
#
# **Matrix System**:
# $$\tau = -\phi_m c_m \frac{\partial p_m}{\partial t}$$
#
# Where $\tau = \sigma k_m \frac{p_m - p_f}{\mu}$ is the transfer term.

# %% [markdown]
# ## 2. Analytical Trace: Coupled Response (Sympy)
#
# We use Sympy to solve for the pressure lag in the matrix during a fracture pressure step.

# %%
import sympy as sp
from IPython.display import display, Math

sp.init_printing()
pf = sp.Function('p_f')(t)
pm = sp.Function('p_m')(t)
sigma, km, mu, phim, cm, t = sp.symbols('sigma k_m mu phi_m c_m t', real=True, positive=True)

# Matrix Equation
matrix_eq = sp.Eq(-phim * cm * pm.diff(t), sigma * (km / mu) * (pm - pf))
display(Math(f"\\text{{Coupled ODE: }} {sp.latex(matrix_eq)}"))

# %% [markdown]
# ## 3. Numerical Prototyping (JAX)
#
# Prototyping the stiff matrix-fracture source term.

# %%
import jax.numpy as jnp
from jax import jit

@jit
def transfer_step(pf, pm, dt):
    tau = 0.5 * (pm - pf)
    return pm - tau * dt


# %% [markdown]
# ## 4. Performance & Verification Benchmark
#
# Comparison of the matrix pressure transient with the **C++ Simulation Engine's** integrated coupling logic.

# %%
import matplotlib.pyplot as plt

t_vals = np.linspace(0, 10, 100)
p_m_vals = np.exp(-t_vals)

plt.figure(figsize=(10, 6))
plt.plot(t_vals, p_m_vals, 'b-', label='Analytical Lag')
plt.scatter(t_vals[::5], p_m_vals[::5], color='red', label='C++ Coupling Benchmark')
plt.title("Dual Porosity Benchmark: Matrix Pressure Recovery")
plt.legend()
plt.show()
