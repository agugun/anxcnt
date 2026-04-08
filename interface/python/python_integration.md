# Python Integration

This directory contains the Python interface for the NumPhys backend, powered by **pybind11**. It bridges the high-performance C++ numerical core with Python\'s interactive ecosystem (Streamlit, Jupyter, etc.).

## Components

- **[`python_bridge.cpp`](python_bridge.cpp)**: The C++ source file that defines the Python bindings. It maps C++ classes (like `Heat1DSimulation`) to Python classes (like `Heat1DImplicit`).
- **[`setup.py`](setup.py)**: The build script used to compile the C++ extension into a shared object (`.so` file on Linux) that Python can import.

## Setup & Build Instructions

To build the Python module, navigate to the `interface/python/` directory and run:

```bash
# Ensure your virtual environment is active
../../.venv/bin/python setup.py build_ext --inplace
```

This will generate a file named `cnt.cpython-*.so` in the current directory.

## Exported API

The following C++ features are exported to Python:

### `Heat1DImplicit` (Class)
- **Constructor**: `Heat1DImplicit(nx, dx, alpha)`
- **`set_initial_condition(ic)`**: Sets the initial temperature field (list of doubles).
- **`set_boundary_conditions(left, right)`**: Sets Dirichlet boundary conditions.
- **`step(dt)`**: Advances the simulation by one time step `dt`.
- **`get_values()`**: Returns the current field values as a list.

## Development Workflow

1.  Modify numerical logic in `../../backend/lib/numerical_methods/` or `../../backend/physics/`.
2.  Update `python_bridge.cpp` if you add new classes or methods.
3.  Re-run the build command.
4.  Verify in the Jupyter Notebook or Streamlit UI.
