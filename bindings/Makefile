# NumPhys Interface Makefile

# Directories
PYTHON_DIR = python
PYTHON_VENV = ../.venv/bin/python

.PHONY: all clean python_inplace

all: python_inplace

python_inplace:
	@echo "Building Python Extension In-place..."
	cd $(PYTHON_DIR) && ../$(PYTHON_VENV) setup.py build_ext --inplace

clean:
	cd $(PYTHON_DIR) && rm -rf build/ *.so
