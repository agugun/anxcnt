import numpy as np
import os
import pyvista as pv
import pandas as pd

def discover_projects(directory="exports"):
    """
    Scans the exports directory for subfolders that follow the Axscnt standard
    (contain a 'vtk' or 'csv' subdirectory).
    """
    projects = []
    if not os.path.exists(directory):
        return projects
        
    for item in os.listdir(directory):
        path = os.path.join(directory, item)
        if os.path.isdir(path):
            # Check if it has a vtk or csv folder
            if os.path.exists(os.path.join(path, 'vtk')) or os.path.exists(os.path.join(path, 'csv')):
                projects.append(item)
    
    return sorted(projects)

def get_available_timesteps(project_name, directory="exports"):
    """
    Returns sorted list of available step numbers for a project by scanning the vtk folder.
    """
    steps = []
    project_path = os.path.join(directory, project_name, "vtk")
    if not os.path.exists(project_path):
        return steps
        
    for f in os.listdir(project_path):
        if f.startswith("step_") and f.endswith(".vti"):
            try:
                # step_0050.vti -> 50
                step_str = f.replace("step_", "").replace(".vti", "")
                steps.append(int(step_str))
            except ValueError:
                continue
    
    return sorted(steps)

def load_sim_snapshot(project_name, step, directory="exports"):
    """
    Loads a specific simulation snapshot using PyVista.
    
    Returns:
        dict: {
            'grid': {'nx', 'ny', 'nz'},
            'spacing': {'dx', 'dy', 'dz'},
            'fields': { name: 3D numpy array }
        }
    """
    filename = f"step_{step:04d}.vti"
    path = os.path.join(directory, project_name, "vtk", filename)
    
    if not os.path.exists(path):
        raise FileNotFoundError(f"VTI file not found: {path}")
        
    mesh = pv.read(path)
    nx, ny, nz = mesh.dimensions
    spacing = mesh.spacing
    
    result = {
        'grid': {'nx': nx, 'ny': ny, 'nz': nz},
        'spacing': {'dx': spacing[0], 'dy': spacing[1], 'dz': spacing[2]},
        'fields': {}
    }
    
    # Extract all point data fields
    for name in mesh.point_data.keys():
        # mesh.point_data[name] is a flat array, reshape to (nz, ny, nx)
        # Note: PyVista/VTK ordering is (X, Y, Z). 
        # mesh.dimensions returns (NX, NY, NZ).
        # We want (NZ, NY, NX) for our internal 2D slicer logic.
        data = mesh.point_data[name]
        result['fields'][name] = data.reshape((nz, ny, nx), order='F') # Fortran order preserves VTK mapping
        
    return result

def get_volumetric_averages(project_name, directory="exports"):
    """
    Attempts to load trend data from results.csv or by scanning snapshots.
    """
    csv_path = os.path.join(directory, project_name, "results.csv")
    if os.path.exists(csv_path):
        try:
            df = pd.read_csv(csv_path)
            # Assuming columns: time, index, x, y, z, field1, field2...
            # We want to group by time/step
            # If the backend results.csv is too detailed (per-cell), we average it here.
            
            # Simple check: does it have a 'time' column?
            if 'time' in df.columns:
                grouped = df.groupby('time').mean()
                time_steps = grouped.index.tolist()
                data_dict = {col: grouped[col].tolist() for col in grouped.columns if col not in ['index', 'x', 'y', 'z']}
                return time_steps, data_dict
        except Exception as e:
            print(f"Error parsing results.csv for {project_name}: {e}")

    # Fallback: Scan snapshots (slower but reliable)
    steps = get_available_timesteps(project_name, directory)
    history = []
    for s in steps:
        try:
            snap = load_sim_snapshot(project_name, s, directory)
            row = {'TimeStep': s}
            for name, arr in snap['fields'].items():
                row[name] = np.mean(arr)
            history.append(row)
        except:
            continue
            
    if not history:
        return [], {}
        
    df = pd.DataFrame(history)
    data_dict = {col: df[col].tolist() for col in df.columns if col != 'TimeStep'}
    return df['TimeStep'].tolist(), data_dict
