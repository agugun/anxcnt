# NumPhys

A numerical physics project incorporating a C++ src and a Python data science/visualization stack.

## Project Structure
- `src/` - C++ kernel and simulation files.
- `bindings/` - Python-C++ bridge and bindings layer.
- `notebook/` - Jupyter notebooks for running experiments and visualizations.
- `web/` - Directory for Streamlit UI components.

## Setup Instructions

Follow these step-by-step instructions to set up the Python environment for the project.

### 1. Prerequisites
- **Python 3.8+** installed on your system.
- **C++ Compiler** (e.g., GCC/Clang or MSVC) for the src.

### 2. Create the Python Virtual Environment
Navigate to the root directory of the project in your terminal, and create a virtual environment:

```bash
python -m venv .venv
```

### 3. Activate the Environment

**On Linux/macOS:**
```bash
source venv/bin/activate
```

**On Windows:**
```bash
.\venv\Scripts\activate
```

### 4. Install Dependencies
Once the environment is activated, install the required packages using `requirements.txt`:

```bash
pip install -r requirements.txt
```

### 5. Start Jupyter (Optional)
To interact with the sample notebooks:

```bash
jupyter notebook
```


### 6. Build the C++ Backend (High Performance)

The src uses a unified build system via a **Makefile**. To compile the Python bridge and the standalone C++ executables:

make all

This will create a `dist/` directory for C++ executables and build the Python bridge:
- `bindings/python/cnt.so`: The Python interface module.
- `dist/heat_1d_implicit`: A standalone C++ executable for the 1D Heat Equation.

To clean the build artifacts, use:
```bash
make clean
```

### 7. Run the Interactive UI (Streamlit)

To launch the real-time simulation dashboard with Plotly visualization:

```bash
streamlit run web/app.py
```
This will open a web interface where you can adjust physical parameters (diffusivity, time step) and observe the thermal evolution instantly.

---

### Building & Debugging Physics Modules
Each physics module (e.g., `heat_1d_implicit`) can be built and debugged dynamically:

1.  **Build a Module**: Press **Ctrl+Shift+B** (or **Run Task** -> **build physics**). You will be prompted to select a module from the list.
2.  **Debug a Module**: Press **F5** (or **Run and Debug** -> **Debug Simulation Physics**). You will be prompted to select a module. VS Code will automatically build the latest version with debug symbols (`-g`) and start the debugger.

3.  **Launch Streamlit UI**: Press **F5** and select **Launch Streamlit UI** from the debug configuration dropdown. This will start the web dashboard using the project's virtual environment and automatically configure the `PYTHONPATH` to find the C++ src.

### Registering New Physics Modules
When you add a new physics folder to `src/physics/` (e.g., `wave_equation`), you must register it in the following files to enable VS Code selection:
- `.vscode/tasks.json`: Add the name to the `options` list of the `selectedCase` input.
- `.vscode/launch.json`: Add the name to the `options` list of the `selectedCase` input.
