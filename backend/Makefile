# NumPhys Backend Makefile

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

heat_sim: $(BUILD_DIR)
	@echo "Building Heat Simulation C++ Executable..."
	$(CXX) $(CXX_FLAGS) physics/heat_1d_implicit/main.cpp -o $(BUILD_DIR)/heat_sim

pressure_sim: $(BUILD_DIR)
	@echo "Building Pressure Simulation C++ Executable..."
	$(CXX) $(CXX_FLAGS) physics/pressure_diffusivity_1d/main.cpp -o $(BUILD_DIR)/pressure_sim

build_physics: $(BUILD_DIR)
	@echo "Building physics module: $(PHYSICS)..."
	$(CXX) $(CXX_FLAGS) -g physics/$(PHYSICS)/main.cpp -o $(BUILD_DIR)/$(PHYSICS)

clean:
	rm -rf $(BUILD_DIR)
