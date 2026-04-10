# NumPhys Consolidated Makefile

# Directories
SRC_DIR = src
BINDINGS_DIR = bindings
BUILD_DIR = dist
EXPORTS_DIR = exports
RESULTS_DIR = results
PYTHON_DIR = $(BINDINGS_DIR)/python
PYTHON_VENV = .venv/bin/python

# Compiler
CXX = g++
CXX_FLAGS = -std=c++14 -O3 -I$(SRC_DIR)

# Targets
.PHONY: all mod bindings clean mod_clean bindings_clean heat_1d_implicit pressure_sim heat_1d_explicit wave heat_2d_implicit heat_3d_implicit mba python_inplace build_physics

all: mod bindings

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)
	mkdir -p $(EXPORTS_DIR) $(RESULTS_DIR)

mod: heat_1d_implicit heat_1d_explicit pressure_sim wave heat_2d_implicit heat_3d_implicit mba

heat_1d_implicit: $(BUILD_DIR)
	@echo "Building Heat Simulation (Implicit) C++ Executable..."
	$(CXX) $(CXX_FLAGS) $(SRC_DIR)/modules/heat_1d_implicit/main.cpp -o $(BUILD_DIR)/heat_1d_implicit

heat_1d_explicit: $(BUILD_DIR)
	@echo "Building Heat Simulation (Explicit) C++ Executable..."
	$(CXX) $(CXX_FLAGS) $(SRC_DIR)/modules/heat_1d_explicit/main.cpp -o $(BUILD_DIR)/heat_1d_explicit

pressure_sim: $(BUILD_DIR)
	@echo "Building Pressure Simulation C++ Executable..."
	$(CXX) $(CXX_FLAGS) $(SRC_DIR)/modules/pressure_diffusivity_1d/main.cpp -o $(BUILD_DIR)/pressure_sim

wave: $(BUILD_DIR)
	@echo "Building Wave Simulation C++ Executable..."
	$(CXX) $(CXX_FLAGS) $(SRC_DIR)/modules/wave_1d/main.cpp -o $(BUILD_DIR)/wave

heat_2d_implicit: $(BUILD_DIR)
	@echo "Building 2D Heat Simulation (Implicit) C++ Executable..."
	$(CXX) $(CXX_FLAGS) $(SRC_DIR)/modules/heat_2d_implicit/main.cpp -o $(BUILD_DIR)/heat_2d_implicit

heat_3d_implicit: $(BUILD_DIR)
	@echo "Building 3D Heat Simulation (Implicit) C++ Executable..."
	$(CXX) $(CXX_FLAGS) $(SRC_DIR)/modules/heat_3d_implicit/main.cpp -o $(BUILD_DIR)/heat_3d_implicit

mba: $(BUILD_DIR)
	@echo "Building Material Balance Simulation C++ Executable..."
	$(CXX) $(CXX_FLAGS) $(SRC_DIR)/modules/material_balance/main.cpp -o $(BUILD_DIR)/mba

build_physics: $(BUILD_DIR)
	@echo "Building physics module: $(PHYSICS)..."
	$(CXX) $(CXX_FLAGS) -g $(SRC_DIR)/modules/$(PHYSICS)/main.cpp -o $(BUILD_DIR)/$(PHYSICS)

bindings: python_inplace

python_inplace:
	@echo "Building Python Extension In-place..."
	cd $(PYTHON_DIR) && ../../$(PYTHON_VENV) setup.py build_ext --inplace

clean: mod_clean bindings_clean

mod_clean:
	@echo "Cleaning C++ Module artifacts..."
	rm -rf $(BUILD_DIR)
	rm -rf $(EXPORTS_DIR) $(RESULTS_DIR)

bindings_clean:
	@echo "Cleaning Python Binding artifacts..."
	cd $(PYTHON_DIR) && rm -rf build/ *.so
