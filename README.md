# AXSCNT 
**Architecture of eXperimental Software for Computational Numerical Techniques**  *Pronounced: /ax-SENT/ (like "Axe-Scent")* is a coding framework designed as a **collaborative bridge** between pure mathematics, domain physics, and software engineering.

*   **Low-Level Transparency**: Prioritize collaboration at the implementation level. Every numerical kernel is structured to be readable for physicists while remaining highly optimized for silicon.
*   **Physics-First Architecture**: Organized by physical disciplines (Thermodynamics, Fluids, Reservoirs), ensuring that the code structure mirrors the mental model of the researcher.
*   **The Analytical Bridge**: Explicitly bridge the gap between "notebook math" and "backend production" through a unified multi-stage workflow.

---

### Workflow
AXSCNT provides a seamless end-to-end pipeline from conceptual derivation to deployable high-performance kernels.

| Stage | Focus | Primary Tools | Role in Ecosystem |
| :--- | :--- | :--- | :--- |
| **1. Formulation** | Theoretical Derivations | **Sympy** | Exact symbolic math & PDE definitions. |
| **2. Prototyping** | Numerical Stability | **JAX / Numpy** | Rapid experimental verification with XLA acceleration. |
| **3. Production** | High Performance | **C++17 / OpenMP** | Deployable, memory-safe, parallel numerical engine. |
| **4. Integration** | Synergy | **axcnt_cpp** | High-speed Python bindings for hybrid workflows. |
| **5. Insight** | Visualization | **Dash / VTK** | Interactive analytics and field monitoring. |

---

### Project Structure

```bash
axscnt/
├── src/
│   ├── lib/            # Foundational Framework & Numerical Utilities (C++)
│   └── modules/        # Domain-Specific Physics Engines (C++)
│       ├── thermodynamics/ 
│       ├── fluids/         
│       ├── reservoir/      
│       ├── wave/           
│       └── pressure/       
├── bindings/           # High-performance C++ connectivity (axcnt_cpp)
├── notebook/           # The Analytical Bridge (Math ➔ Python ➔ C++)
│   ├── thermodynamics/ 
│   ├── fluids/         
│   ├── reservoir/      
│   └── ...
├── web/                # Collaborative Visual Analytics (Python/Dash)
├── benchmark/          # Performance & Scalability Benchmarking
├── dist/               # Production-ready C++ Executables
└── exports/            # Standardized Simulation Outputs (VTK, JSON, CSV)
```

---

### Internal Architecture

The C++ core is strictly segregated by namespaces to ensure mathematical consistency and infinite modularity.

#### The Numerical Engine (`top` & `num`)
- **Architectural Contracts (`top`)**: Defines the abstract interfaces (`IState`, `IModel`, `ISolver`) that enforce the AXSCNT "Simulation Contract".
- **Numerical Core (`num`)**: Contains vectorized linear solvers (BiCGSTAB, Tridiagonal), non-linear Newton-Raphson schemes, and sparse matrix implementations.
- **Infrastructure (`utl`)**: Handles cross-cutting concerns like hierarchical configuration, standardized logging, and filesystem orchestration.

#### Physics Module Anatomy (`mod`)
Every physical module implements the framework contract through four key components:

1.  **`model.hpp`**: Defines the physical constants and the **Governing Equation** discretization.
2.  **`state.hpp`**: A structured snapshot of the field data at any point in time.
3.  **`simulation.hpp`**: A factory/builder that assembles the complex object graph (Grid ➔ Discretizer ➔ Linearizer ➔ Solver).
4.  **`main.cpp`**: The execution orchestrator for standalone deployment.

---

### Quick Start

#### 1. Setup the Analytical Environment
```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

#### 2. Compile the Production Engine
```bash
make all
```

#### 3. Run Your First Simulation
```bash
# Run a 2D Reservoir Simulation
./dist/reservoir_2d input/reservoir_2d.txt

# Launch the Visual Analytics Dashboard
python web/app.py
```

---

### Openness
AXSCNT is built from & for the community.
- **Copy-Paste Encouraged**: Take any kernel, discretization strategy, or utility for your own projects.
- **No Restrictions**: Use it for research, production, or education without friction.
- **Contributions**: We welcome domain experts to implement new physics modules using our standardized `top` interfaces.

---
