import sys
import os
import unittest

# Add bindings path to sys.path
binding_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'bindings', 'python'))
sys.path.append(binding_path)

try:
    import axcnt_cpp
except ImportError as e:
    print(f"Failed to import axcnt_cpp: {e}")
    sys.exit(1)

class TestPythonBindings(unittest.TestCase):
    def test_import(self):
        self.assertIsNotNone(axcnt_cpp, "axcnt_cpp module should be loaded")

    def test_heat_1d(self):
        # Test basic initialization of Heat1D model
        nx = 100
        # Create dummy conductance and storage arrays
        T_cond = [0.05] * (nx - 1)
        storage = [1.0] * nx
        model = axcnt_cpp.Heat1D(T_cond, storage, 0.0, 0.0)
        self.assertIsNotNone(model)
        
        # Test setting initial condition
        u0 = [0.0] * nx
        u0[50] = 1.0
        model.set_initial_condition(u0)
        
        # Test taking a step
        model.step(0.01)
        
        values = model.get_values()
        self.assertEqual(len(values), nx)
        self.assertAlmostEqual(values[0], 0.0)
        self.assertAlmostEqual(values[-1], 0.0)

if __name__ == '__main__':
    unittest.main()
