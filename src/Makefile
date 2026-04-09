# NumPhys Source Makefile

# Directories
BUILD_DIR = build

# Compiler
CXX = g++
CXX_FLAGS = -std=c++14 -O3 -I.

# Targets
.PHONY: all clean heat_sim pressure_sim

all: heat_sim pressure_sim

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)
	mkdir -p exports

heat_sim: $(BUILD_DIR)
	@echo "Building Heat Simulation C++ Executable..."
	$(CXX) $(CXX_FLAGS) modules/heat_1d_implicit/main.cpp -o $(BUILD_DIR)/heat_sim

heat: $(BUILD_DIR)
	@echo "Building Heat Simulation C++ Executable..."
	$(CXX) $(CXX_FLAGS) modules/heat/main.cpp -o $(BUILD_DIR)/heat

pressure_sim: $(BUILD_DIR)
	@echo "Building Pressure Simulation C++ Executable..."
	$(CXX) $(CXX_FLAGS) modules/pressure_diffusivity_1d/main.cpp -o $(BUILD_DIR)/pressure_sim

wave: $(BUILD_DIR)
	@echo "Building Wave Simulation C++ Executable..."
	$(CXX) $(CXX_FLAGS) modules/wave_1d/main.cpp -o $(BUILD_DIR)/wave

heat_2d: $(BUILD_DIR)
	@echo "Building 2D Heat Simulation C++ Executable..."
	$(CXX) $(CXX_FLAGS) modules/heat_2d_implicit/main.cpp -o $(BUILD_DIR)/heat_2d

mba: $(BUILD_DIR)
	@echo "Building Material Balance Simulation C++ Executable..."
	$(CXX) $(CXX_FLAGS) modules/material_balance/main.cpp -o $(BUILD_DIR)/mba

build_physics: $(BUILD_DIR)
	@echo "Building physics module: $(PHYSICS)..."
	$(CXX) $(CXX_FLAGS) -g modules/$(PHYSICS)/main.cpp -o $(BUILD_DIR)/$(PHYSICS)

clean:
	rm -rf $(BUILD_DIR)
