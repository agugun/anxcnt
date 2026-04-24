# NumPhys Consolidated Makefile

# Directories
SRC_DIR = src
BINDINGS_DIR = bindings
BUILD_DIR = dist
EXPORTS_DIR = exports
PYTHON_DIR = $(BINDINGS_DIR)/python
PYTHON_VENV = .venv/bin/python

# Compiler
CXX = g++
CXX_FLAGS = -std=c++17 -O3 -fopenmp -I$(SRC_DIR)

# Targets
.PHONY: all mod bindings clean mod_clean bindings_clean heat_1d_implicit pressure_sim heat_1d_explicit wave_1d heat_2d_implicit heat_3d_implicit mba reservoir_1d reservoir_2d reservoir_3d reservoir_dual_2d reservoir_oil_gas_2d reservoir_black_oil_2d reservoir_black_oil_3d burgers_fdm fluid_fem wave_2d python_inplace build_physics

all: mod bindings

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

mod: heat_1d_implicit heat_1d_explicit pressure_sim wave_1d heat_2d_implicit heat_3d_implicit mba reservoir_1d reservoir_2d reservoir_3d reservoir_dual_2d reservoir_oil_gas_2d reservoir_black_oil_2d reservoir_black_oil_3d burgers_fdm fluid_fem wave_2d oscillator

heat_1d_implicit: $(BUILD_DIR)
	@echo "Building Heat Simulation (Implicit) C++ Executable..."
	$(CXX) $(CXX_FLAGS) $(SRC_DIR)/modules/thermodynamics/heat/1d_implicit/main.cpp -o $(BUILD_DIR)/heat_1d_implicit

heat_1d_explicit: $(BUILD_DIR)
	@echo "Building Heat Simulation (Explicit) C++ Executable..."
	$(CXX) $(CXX_FLAGS) $(SRC_DIR)/modules/thermodynamics/heat/1d_explicit/main.cpp -o $(BUILD_DIR)/heat_1d_explicit

pressure_sim: $(BUILD_DIR)
	@echo "Building Pressure Simulation C++ Executable..."
	$(CXX) $(CXX_FLAGS) $(SRC_DIR)/modules/pressure/1d/main.cpp -o $(BUILD_DIR)/pressure_sim

wave_1d: $(BUILD_DIR)
	@echo "Building 1D Wave Simulation C++ Executable..."
	$(CXX) $(CXX_FLAGS) $(SRC_DIR)/modules/wave/1d/main.cpp -o $(BUILD_DIR)/wave_1d

wave_2d: $(BUILD_DIR)
	@echo "Building 2D Wave Equation Simulation (FDM) C++ Executable..."
	$(CXX) $(CXX_FLAGS) $(SRC_DIR)/modules/wave/2d/main.cpp -o $(BUILD_DIR)/wave_2d

heat_2d_implicit: $(BUILD_DIR)
	@echo "Building 2D Heat Simulation (Implicit) C++ Executable..."
	$(CXX) $(CXX_FLAGS) $(SRC_DIR)/modules/thermodynamics/heat/2d_implicit/main.cpp -o $(BUILD_DIR)/heat_2d_implicit

heat_3d_implicit: $(BUILD_DIR)
	@echo "Building 3D Heat Simulation (Implicit) C++ Executable..."
	$(CXX) $(CXX_FLAGS) $(SRC_DIR)/modules/thermodynamics/heat/3d_implicit/main.cpp -o $(BUILD_DIR)/heat_3d_implicit

mba: $(BUILD_DIR)
	@echo "Building Material Balance Simulation C++ Executable..."
	$(CXX) $(CXX_FLAGS) $(SRC_DIR)/modules/reservoir/mba/main.cpp -o $(BUILD_DIR)/mba

reservoir_1d: $(BUILD_DIR)
	@echo "Building 1D Reservoir Simulation C++ Executable..."
	$(CXX) $(CXX_FLAGS) $(SRC_DIR)/modules/reservoir/1d/main.cpp -o $(BUILD_DIR)/reservoir_1d

reservoir_2d: $(BUILD_DIR)
	@echo "Building 2D Reservoir Simulation C++ Executable..."
	$(CXX) $(CXX_FLAGS) $(SRC_DIR)/modules/reservoir/2d/main.cpp -o $(BUILD_DIR)/reservoir_2d

reservoir_3d: $(BUILD_DIR)
	@echo "Building 3D Reservoir Simulation C++ Executable..."
	$(CXX) $(CXX_FLAGS) $(SRC_DIR)/modules/reservoir/3d/main.cpp -o $(BUILD_DIR)/reservoir_3d

reservoir_dual_2d: $(BUILD_DIR)
	@echo "Building 2D Dual-Phase Reservoir Simulation C++ Executable..."
	$(CXX) $(CXX_FLAGS) $(SRC_DIR)/modules/reservoir/dual_2d/main.cpp -o $(BUILD_DIR)/reservoir_dual_2d

reservoir_oil_gas_2d: $(BUILD_DIR)
	@echo "Building 2D Oil-Gas Reservoir Simulation C++ Executable..."
	$(CXX) $(CXX_FLAGS) $(SRC_DIR)/modules/reservoir/oil_gas_2d/main.cpp -o $(BUILD_DIR)/reservoir_oil_gas_2d

reservoir_black_oil_2d: $(BUILD_DIR)
	@echo "Building 2D Black Oil Reservoir Simulation C++ Executable..."
	$(CXX) $(CXX_FLAGS) $(SRC_DIR)/modules/reservoir/black_oil_2d/main.cpp -o $(BUILD_DIR)/reservoir_black_oil_2d

reservoir_black_oil_3d: $(BUILD_DIR)
	@echo "Building 3D Black Oil Reservoir Simulation C++ Executable..."
	$(CXX) $(CXX_FLAGS) $(SRC_DIR)/modules/reservoir/black_oil_3d/main.cpp -o $(BUILD_DIR)/reservoir_black_oil_3d

burgers_fdm: $(BUILD_DIR)
	@echo "Building 1D Burgers' Equation Simulation (FDM) C++ Executable..."
	$(CXX) $(CXX_FLAGS) $(SRC_DIR)/modules/fluids/burgers/main.cpp -o $(BUILD_DIR)/burgers_fdm

fluid_fem: $(BUILD_DIR)
	@echo "Building 2D Fluid Dynamics Simulation (FEM) C++ Executable..."
	$(CXX) $(CXX_FLAGS) $(SRC_DIR)/modules/fluids/fluid_dynamics/main.cpp -o $(BUILD_DIR)/fluid_fem

oscillator: $(BUILD_DIR)
	@echo "Building Harmonic Oscillator Pilot Case C++ Executable..."
	$(CXX) $(CXX_FLAGS) $(SRC_DIR)/modules/oscillator/main.cpp -o $(BUILD_DIR)/oscillator


build_physics: $(BUILD_DIR)
	@echo "Building physics module: $(PHYSICS)..."
	$(CXX) $(CXX_FLAGS) -g $(SRC_DIR)/modules/$(PHYSICS)/main.cpp -o $(BUILD_DIR)/$(notdir $(PHYSICS))

bindings: python_inplace

python_inplace:
	@echo "Building Python Extension In-place..."
	cd $(PYTHON_DIR) && ../../$(PYTHON_VENV) setup.py build_ext --inplace

build_test_math: $(BUILD_DIR)
	@echo "Building Math Verification Test..."
	$(CXX) $(CXX_FLAGS) tests/verify_math.cpp -o $(BUILD_DIR)/verify_math

build_test_solvers: $(BUILD_DIR)
	@echo "Building Solver Verification Test..."
	$(CXX) $(CXX_FLAGS) tests/verify_solvers.cpp -o $(BUILD_DIR)/verify_solvers

build_test_integrators: $(BUILD_DIR)
	@echo "Building Integrator Verification Test..."
	$(CXX) $(CXX_FLAGS) tests/verify_integrators.cpp -o $(BUILD_DIR)/verify_integrators

test_math: build_test_math
	@echo "Running Math Verification Tests..."
	./$(BUILD_DIR)/verify_math

test_solvers: build_test_solvers
	@echo "Running Solver Verification Tests..."
	./$(BUILD_DIR)/verify_solvers

test_integrators: build_test_integrators
	@echo "Running Integrator Verification Tests..."
	./$(BUILD_DIR)/verify_integrators

test_python: bindings
	@echo "Running Python Binding Tests..."
	$(PYTHON_VENV) tests/test_python_bindings.py

test: test_math test_solvers test_integrators test_python

build_benchmark: $(BUILD_DIR)
	@echo "Building Solver Performance Benchmark..."
	$(CXX) $(CXX_FLAGS) benchmark/benchmark_solvers.cpp -o $(BUILD_DIR)/benchmark_solvers

benchmark: build_benchmark
	@echo "Running Solver Performance Benchmarks..."
	./$(BUILD_DIR)/benchmark_solvers

clean: mod_clean bindings_clean

mod_clean:
	@echo "Cleaning C++ Module artifacts..."
	rm -rf $(BUILD_DIR)
	rm -rf $(EXPORTS_DIR)

bindings_clean:
	@echo "Cleaning Python Binding artifacts..."
	cd $(PYTHON_DIR) && rm -rf build/ *.so
