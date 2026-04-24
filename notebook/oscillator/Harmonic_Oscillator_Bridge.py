# %% [markdown]
# # Harmonic Oscillator

# %% [markdown]
# ### Analytical Formulation
# 
# The harmonic oscillator represents a fundamental physical system where a mass $m$ is subject to a restoring force proportional to its displacement $x$ (Hooke's Law: $F = -kx$, where $k$ is the spring constant). This system models a wide variety of phenomena, from simple pendulums to molecular vibrations.
#
# By applying Newton's Second Law ($F = ma$), we define the undamped harmonic oscillator equation:
# $$ m \frac{d^2x}{dt^2} + k x = 0 $$
# With initial conditions $x(0) = 1$, $x'(0) = 0$.

# %%
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

# Fallback for display() if running as a raw Python script outside of Jupyter/IPython
try:
    from IPython.display import display
except ImportError:
    def display(x): print(x)

# Shared Parameters for testing later
M = 1.0
K = 4.0
T_MAX = 10.0
DT = 0.05

# 1. Define symbols
t, m, k = sp.symbols('t m k', positive=True)
x = sp.Function('x')(t)

# 2. Define the differential equation
ode = sp.Eq(m * x.diff(t, 2) + k * x, 0)
display(ode)

# 3. Solve with initial conditions: x(0) = 1, x'(0) = 0
ics = {x.subs(t, 0): 1, x.diff(t).subs(t, 0): 0}
analytical_sol = sp.dsolve(ode, x, ics=ics)
display(analytical_sol)

# 4. Create a fast, callable numerical function from the SymPy expression
x_exact_func = sp.lambdify((t, m, k), analytical_sol.rhs, "numpy")

# %% [markdown]
# ### Numerical Implementation
# To solve this numerically, we translate the 2nd-order ODE into a system of 1st-order ODEs.
# We also implement a non-linear Duffing oscillator to demonstrate breaking analytical limits.

# %%
def solve_numerical(m_val, k_val, t_max, dt, method='rk4'):
    """Numerical solver for the linear harmonic oscillator."""
    times = np.arange(0, t_max, dt)
    n_steps = len(times)
    
    x_num = np.zeros(n_steps)
    v_num = np.zeros(n_steps)
    
    x_num[0] = 1.0  # Initial position
    v_num[0] = 0.0  # Initial velocity
    
    for i in range(1, n_steps):
        x_prev = x_num[i-1]
        v_prev = v_num[i-1]
        
        if method == 'euler':
            # Explicit Euler (Usually unstable for oscillatory systems)
            x_num[i] = x_prev + dt * v_prev
            v_num[i] = v_prev - dt * (k_val / m_val * x_prev)
            
        elif method == 'rk4':
            # Runge-Kutta 4 (Highly stable, 4th order accuracy)
            k1_x = v_prev
            k1_v = -(k_val / m_val) * x_prev
            
            k2_x = v_prev + 0.5 * dt * k1_v
            k2_v = -(k_val / m_val) * (x_prev + 0.5 * dt * k1_x)
            
            k3_x = v_prev + 0.5 * dt * k2_v
            k3_v = -(k_val / m_val) * (x_prev + 0.5 * dt * k2_x)
            
            k4_x = v_prev + dt * k3_v
            k4_v = -(k_val / m_val) * (x_prev + dt * k3_x)
            
            x_num[i] = x_prev + (dt / 6.0) * (k1_x + 2*k2_x + 2*k3_x + k4_x)
            v_num[i] = v_prev + (dt / 6.0) * (k1_v + 2*k2_v + 2*k3_v + k4_v)

        elif method == 'implicit_euler':
            # Implicit Euler (Unconditionally stable, but highly dissipative)
            v_num[i] = (v_prev - dt * (k_val / m_val) * x_prev) / (1.0 + dt**2 * (k_val / m_val))
            x_num[i] = x_prev + dt * v_num[i]
            
    return times, x_num

def solve_duffing(m_val, k_val, c, alpha, t_max, dt):
    """Numerical solver for the NON-LINEAR Duffing oscillator."""
    # Solves dx/dt = v; dv/dt = -(k/m) x - (c/m) v - (alpha/m) x^3
    times = np.arange(0, t_max, dt)
    n_steps = len(times)
    
    x_num = np.zeros(n_steps)
    v_num = np.zeros(n_steps)
    
    x_num[0] = 1.0  # Initial position
    v_num[0] = 0.0  # Initial velocity
    
    for i in range(1, n_steps):
        x_prev = x_num[i-1]
        v_prev = v_num[i-1]
        
        def dvdt(x_val, v_val):
            # The complex physics that SymPy can't easily solve
            return - (k_val / m_val) * x_val - (c / m_val) * v_val - (alpha / m_val) * (x_val**3)
            
        k1_x = v_prev
        k1_v = dvdt(x_prev, v_prev)
        
        k2_x = v_prev + 0.5 * dt * k1_v
        k2_v = dvdt(x_prev + 0.5 * dt * k1_x, v_prev + 0.5 * dt * k1_v)
        
        k3_x = v_prev + 0.5 * dt * k2_v
        k3_v = dvdt(x_prev + 0.5 * dt * k2_x, v_prev + 0.5 * dt * k2_v)
        
        k4_x = v_prev + dt * k3_v
        k4_v = dvdt(x_prev + dt * k3_x, v_prev + dt * k3_v)
        
        x_num[i] = x_prev + (dt / 6.0) * (k1_x + 2*k2_x + 2*k3_x + k4_x)
        v_num[i] = v_prev + (dt / 6.0) * (k1_v + 2*k2_v + 2*k3_v + k4_v)
            
    return times, x_num

# %% [markdown]
# ### Testing & Verification
# We plot the Analytical ground truth against the Numerical Python prototypes, and against the production C++ engine results.

# %%
# Get Analytical Solution
t_exact = np.linspace(0, T_MAX, 500)
x_exact = x_exact_func(t_exact, M, K)

# Get Numerical Solutions
t_euler, x_euler = solve_numerical(M, K, T_MAX, DT, method='euler')
t_implicit, x_implicit = solve_numerical(M, K, T_MAX, DT, method='implicit_euler')
t_rk4, x_rk4 = solve_numerical(M, K, T_MAX, DT, method='rk4')
t_duffing, x_duffing = solve_duffing(M, K, c=0.2, alpha=1.0, t_max=T_MAX, dt=DT)

# Plot 1: Python Prototype Verification
plt.figure(figsize=(12, 6))
plt.plot(t_exact, x_exact, 'k-', linewidth=3, alpha=0.5, label='Analytical (SymPy Ground Truth)')
plt.plot(t_euler, x_euler, 'r--', label='Numerical (Explicit Euler - Fails)')
plt.plot(t_implicit, x_implicit, 'g--', label='Numerical (Implicit Euler - Dissipates)')
plt.plot(t_rk4, x_rk4, 'b.', label='Numerical (RK4 - Stable)')
plt.title('Validation: Harmonic Oscillator Analytical vs Numerical')
plt.xlabel('Time (t)')
plt.ylabel('Position (x)')
plt.legend()
plt.grid(True)
plt.show()

# Error Analysis 
x_exact_at_rk4_times = x_exact_func(t_rk4, M, K)
l2_error = np.sqrt(np.mean((x_rk4 - x_exact_at_rk4_times)**2))
print(f"Mathematical Validation -> RK4 L2 Error (dt={DT}): {l2_error:.2e}")

# Plot 2: Non-Linear Limits
plt.figure(figsize=(12, 6))
plt.plot(t_exact, x_exact, 'k--', alpha=0.3, linewidth=2, label='Linear Ground Truth (Analytical)')
plt.plot(t_duffing, x_duffing, 'g-', linewidth=2, label='Non-Linear Duffing Oscillator (Numerical RK4)')
plt.title('Breaking the Limits: Adding Non-Linearity & Damping')
plt.xlabel('Time (t)')
plt.ylabel('Position (x)')
plt.legend()
plt.grid(True)
plt.show()

# Plot 3: Golden Prototype (C++ Execution)
cpp_output_path = "../../exports/oscillator_output.csv"

if os.path.exists(cpp_output_path):
    df_cpp = pd.read_csv(cpp_output_path)
    
    plt.figure(figsize=(12, 6))
    plt.plot(t_exact, x_exact, 'k-', linewidth=4, alpha=0.3, label='Analytical (SymPy Ground Truth)')
    plt.plot(t_rk4, x_rk4, 'b--', linewidth=2, label='Python Prototype (RK4)')
    plt.plot(df_cpp['Time'], df_cpp['Position'], 'ro', markersize=4, label='C++ Production Engine (Implicit Euler)')
    
    plt.title('Closing the Loop: C++ Execution vs Analytical Ground Truth')
    plt.xlabel('Time (t)')
    plt.ylabel('Position (x)')
    plt.legend()
    plt.grid(True)
    plt.show()
    
    # Calculate Error between C++ and Analytical
    x_exact_cpp_times = x_exact_func(df_cpp['Time'].values, M, K)
    cpp_error = np.sqrt(np.mean((df_cpp['Position'].values - x_exact_cpp_times)**2))
    print(f"C++ Engine L2 Error (dt={DT}): {cpp_error:.2e}")
else:
    print("Run `make oscillator` and `./dist/oscillator` to generate the C++ output first!")
