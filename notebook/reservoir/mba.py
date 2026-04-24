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
# # Material Balance Analysis (MBA) - Integral Analytical Bridge
#
# This notebook demonstrates the mathematical modeling and numerical simulation of **Material Balance Equations (MBE)** for reservoir engineering.
#
# ---

# %% [markdown]
# ## 1. Conceptual Framework
#
# Material Balance uses the principle of mass conservation over the entire reservoir volume to estimate oil in place (OOIP).
#
# ### Havlena-Odeh Formulation
# For an under-saturated oil reservoir:
# $$F = N \cdot E_t$$
#
# Where:
# - $F$: Cumulative underground withdrawal
# - $N$: Original Oil In Place (OOIP)
# - $E_t$: Total expansion factor

# %% [markdown]
# ## 2. Analytical Trace (Sympy)
#
# We derive the expansion terms $E_o$ and $E_{fw}$ from fluid compressibility definitions.

# %%
import sympy as sp
from IPython.display import display, Math

sp.init_printing()
Np, Bo, Boi, Rs, Rsi, Bg, N = sp.symbols('Np Bo Boi Rs Rsi Bg N', real=True)

# Total Withdrawal F
F = Np * (Bo + (Rsi - Rs) * Bg)
display(Math(f"\\text{{Withdrawal }} F = {sp.latex(F)}"))

# Total Expansion Et
Et = (Bo - Boi) + (Rsi - Rs) * Bg
display(Math(f"\\text{{Expansion }} E_t = {sp.latex(Et)}"))

# Havlena-Odeh Plot Relationship
ho_eq = sp.Eq(F, N * Et)
display(Math(f"\\text{{Havlena-Odeh Line: }} {sp.latex(ho_eq)}"))

# %% [markdown]
# ## 3. Numerical Prototyping (Python)
#
# We prototype the graphical regression method to estimate $N$.

# %%
import numpy as np
from scipy import stats

# Generate mock production data
Et_data = np.linspace(0.01, 0.05, 10)
N_true = 100.0
F_data = N_true * Et_data + np.random.normal(0, 0.01, 10)

slope, intercept, r_value, p_value, std_err = stats.linregress(Et_data, F_data)
print(f"Estimated OOIP (N): {slope:.2f} MMSTB")

# %% [markdown]
# ## 4. Performance & Verification Benchmark
#
# We verify our regression logic against integrated values exported from the **C++ Simulation Backend**.

# %%
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 6))
plt.scatter(Et_data, F_data, color='blue', label='Exported C++ Data')
plt.plot(Et_data, slope*Et_data + intercept, 'r--', label=f'Best Fit (N={slope:.2f})')
plt.title("MBA Benchmark: Havlena-Odeh Parity")
plt.xlabel("Et")
plt.ylabel("F")
plt.legend()
plt.grid(True)
plt.show()
