# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Axscnt Reservoir Simulation: 3D Visualization Framework
#
# This notebook demonstrates how to use the standardized **Binary + JSON** data format to load simulation results and generate interactive 3D plots in Python.
#
# ### Prerequisites
# Ensure you have the following packages installed in your `.venv`:
# ```bash
# pip install numpy matplotlib pyvista trame trame-vtk trame-vuetify nest_asyncio2
# ```

# %%
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# 1. Setup paths to import our framework utilities
sys.path.append(os.path.abspath("../notebook/lib"))
from data_loader import load_sim_snapshot, get_available_timesteps

print("Framework loader initialized.")

# %% [markdown]
# ## 1. Load Simulation Results
# The simulation exports snapshots in high-performance Binary + JSON format.

# %%
try:
    steps = get_available_timesteps("../exports")
    if not steps:
        raise FileNotFoundError("No snapshots found in exports/ directory.")
    
    print(f"Discovered {len(steps)} timesteps: {steps[:5]} ... {steps[-1]}")
    
    # Load initial step
    data = load_sim_snapshot(f"../exports/res_t{steps[0]}")
    print(f"Successfully loaded grid: {data['grid']}")
    print(f"Available fields: {list(data['fields'].keys())}")
except Exception as e:
    print(f"Setup Error: {e}")
    print("\nTIP: Run simulation with --json --csv --vtk flags first.")

# %% [markdown]
# ## 2. Statistical Analysis
# Analyze field-wide trends for the selected snapshot.

# %%
if 'data' in locals():
    pressure = data['fields']['Pressure']
    
    print(f"Max Pressure: {np.max(pressure):.2f} psi")
    print(f"Min Pressure: {np.min(pressure):.2f} psi")
    
    plt.figure(figsize=(8, 4))
    plt.hist(pressure.flatten(), bins=50, color='skyblue', edgecolor='black')
    plt.title("Pressure Distribution")
    plt.xlabel("Pressure (psi)")
    plt.show()

# %% [markdown]
# ## 3. 2D Slice View
# A standard Matplotlib view for quick profiling.

# %%
if 'data' in locals():
    layer = 0 # Top layer
    plt.figure(figsize=(10, 8))
    plt.imshow(pressure[layer, :, :], cmap='viridis', origin='lower')
    plt.colorbar(label='Pressure (psi)')
    plt.title(f"Reservoir Pressure Profile - Layer {layer}")
    plt.show()

# %% [markdown]
# ## 4. Interactive 3D Visualization with Time Slider
# Watch the reservoir evolve in 3D using PyVista and Trame.

# %%
try:
    import pyvista as pv
    
    # 1. Initialize mesh
    mesh = pv.ImageData(
        dimensions=(data['grid']['nx'], data['grid']['ny'], data['grid']['nz']),
        spacing=(data['spacing']['dx'], data['spacing']['dy'], data['spacing']['dz'])
    )
    mesh.point_data["Pressure"] = data['fields']['Pressure'].flatten()
    
    plotter = pv.Plotter(notebook=True)
    plotter.add_mesh(mesh, scalars="Pressure", cmap="viridis", 
                     clim=[np.min(mesh["Pressure"]), np.max(mesh["Pressure"])])
    
    # 2. Callback for Time Slider
    def update_snapshot(value):
        step = int(value)
        try:
            new_data = load_sim_snapshot(f"../exports/res_t{step}")
            mesh.point_data["Pressure"] = new_data['fields']['Pressure'].flatten()
            plotter.render()
        except Exception as e:
            print(f"Error updating snapshot for step {step}: {e}")

    plotter.add_slider_widget(
        callback=update_snapshot, 
        rng=[min(steps), max(steps)], 
        value=min(steps), 
        title="Simulation Step",
        style='modern'
    )
    
    plotter.show()
except Exception as e:
    print(f"Visualization error: {e}")
