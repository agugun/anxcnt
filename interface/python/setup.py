from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext
import os

# Define the C++ extension module
ext_module = Pybind11Extension(
    "cnt",
    ["python_bridge.cpp"],
    include_dirs=[
        os.path.abspath("../../backend"), # To find case/ and lib/
    ],
    cxx_std=14,
)

setup(
    name="cnt",
    version="0.1.0",
    ext_modules=[ext_module],
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
)
