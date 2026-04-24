import nbformat as nbf

nb = nbf.v4.new_notebook()

nb.cells.append(nbf.v4.new_markdown_cell("""\
# The Analytical Bridge: Harmonic Oscillator
This notebook demonstrates the transition from a pure analytical formulation to a robust numerical solver using the simple Harmonic Oscillator.
"""))

nb.cells.append(nbf.v4.new_markdown_cell("""\
## Phase 1: The Analytical Ground Truth (SymPy)
We define the undamped harmonic oscillator: $\\frac{d^2x}{dt^2} + \omega^2 x = 0$
"""))

nb.cells.append(nbf.v4.new_code_cell("""\
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

# 1. Define symbols
t, omega = sp.symbols('t omega', positive=True)
x = sp.Function('x')(t)

# 2. Define the differential equation
ode = sp.Eq(x.diff(t, 2) + omega**2 * x, 0)
display(ode)

# 3. Solve with initial conditions: x(0) = 1, x'(0) = 0
ics = {x.subs(t, 0): 1, x.diff(t).subs(t, 0): 0}
analytical_sol = sp.dsolve(ode, x, ics=ics)
display(analytical_sol)

# 4. Create a callable numerical function from the SymPy expression
x_exact_func = sp.lambdify((t, omega), analytical_sol.rhs, "numpy")
"""))

nb.cells.append(nbf.v4.new_markdown_cell("""\
## Phase 2: The Numerical Prototype (Explicit Euler & RK4)
We translate the 2nd-order ODE into a system of 1st-order ODEs:
1. $\\frac{dx}{dt} = v$
2. $\\frac{dv}{dt} = -\omega^2 x$
"""))

nb.cells.append(nbf.v4.new_code_cell("""\
def solve_numerical(omega_val, t_max, dt, method='rk4'):
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
            # Explicit Euler
            x_num[i] = x_prev + dt * v_prev
            v_num[i] = v_prev - dt * (omega_val**2 * x_prev)
            
        elif method == 'rk4':
            # Runge-Kutta 4
            k1_x = v_prev
            k1_v = -(omega_val**2) * x_prev
            
            k2_x = v_prev + 0.5 * dt * k1_v
            k2_v = -(omega_val**2) * (x_prev + 0.5 * dt * k1_x)
            
            k3_x = v_prev + 0.5 * dt * k2_v
            k3_v = -(omega_val**2) * (x_prev + 0.5 * dt * k2_x)
            
            k4_x = v_prev + dt * k3_v
            k4_v = -(omega_val**2) * (x_prev + dt * k3_x)
            
            x_num[i] = x_prev + (dt / 6.0) * (k1_x + 2*k2_x + 2*k3_x + k4_x)
            v_num[i] = v_prev + (dt / 6.0) * (k1_v + 2*k2_v + 2*k3_v + k4_v)
            
    return times, x_num
"""))

nb.cells.append(nbf.v4.new_markdown_cell("""\
## Phase 3: The Bridge & Validation
We plot the Analytical ground truth against the Numerical prototype. We will observe that Explicit Euler fails (gains energy) for this oscillatory system, while RK4 is stable.
"""))

nb.cells.append(nbf.v4.new_code_cell("""\
# Shared Parameters
OMEGA = 2.0
T_MAX = 10.0
DT = 0.05

# Get Analytical Solution
t_exact = np.linspace(0, T_MAX, 500)
x_exact = x_exact_func(t_exact, OMEGA)

# Get Numerical Solutions
t_euler, x_euler = solve_numerical(OMEGA, T_MAX, DT, method='euler')
t_rk4, x_rk4 = solve_numerical(OMEGA, T_MAX, DT, method='rk4')

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(t_exact, x_exact, 'k-', linewidth=2, label='Analytical (SymPy)')
plt.plot(t_euler, x_euler, 'r--', label='Numerical (Explicit Euler)')
plt.plot(t_rk4, x_rk4, 'b.', label='Numerical (RK4)')

plt.title('Harmonic Oscillator: Analytical vs Numerical')
plt.xlabel('Time (t)')
plt.ylabel('Position (x)')
plt.legend()
plt.grid(True)
plt.show()

# Error Analysis
x_exact_at_rk4_times = x_exact_func(t_rk4, OMEGA)
l2_error = np.sqrt(np.mean((x_rk4 - x_exact_at_rk4_times)**2))
print(f"RK4 L2 Error (dt={DT}): {l2_error:.2e}")
"""))

nb.cells.append(nbf.v4.new_markdown_cell("""\
## Phase 4: Breaking the Limits (Non-linear Duffing Oscillator)
Now that we trust our RK4 solver, we can easily introduce a non-linear term (like a hardening spring: $x^3$) or damping ($c \cdot v$). The analytical approach with SymPy fails to give a simple closed-form solution here, but our numerical engine handles it effortlessly.
"""))

nb.cells.append(nbf.v4.new_code_cell("""\
def solve_duffing(omega_val, c, alpha, t_max, dt):
    # Solves dx/dt = v; dv/dt = -omega^2 x - c v - alpha x^3
    times = np.arange(0, t_max, dt)
    n_steps = len(times)
    
    x_num = np.zeros(n_steps)
    v_num = np.zeros(n_steps)
    
    x_num[0] = 1.0  # Initial position
    v_num[0] = 0.0  # Initial velocity
    
    for i in range(1, n_steps):
        x_prev = x_num[i-1]
        v_prev = v_num[i-1]
        
        # RK4 Implementation for Duffing
        def dvdt(x_val, v_val):
            return - (omega_val**2) * x_val - c * v_val - alpha * (x_val**3)
            
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

t_duffing, x_duffing = solve_duffing(OMEGA, c=0.2, alpha=1.0, t_max=T_MAX, dt=DT)

plt.figure(figsize=(10, 6))
plt.plot(t_exact, x_exact, 'k--', alpha=0.5, label='Linear (Analytical)')
plt.plot(t_duffing, x_duffing, 'g-', linewidth=2, label='Non-Linear Duffing (Numerical RK4)')
plt.title('Breaking the Limits: Adding Non-Linearity & Damping')
plt.xlabel('Time (t)')
plt.ylabel('Position (x)')
plt.legend()
plt.grid(True)
plt.show()
"""))

with open('/home/gugun/repo/axscnt/notebook/Harmonic_Oscillator_Bridge.ipynb', 'w') as f:
    nbf.write(nb, f)

print("Notebook generated successfully.")
