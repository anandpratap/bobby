import sys, os
import numpy as np
import unittest
import numpy.testing as npt

sys.path.append(os.path.join(os.path.dirname(__file__), "../"))
from src.shapefunctions import ShapeFunction1D
from src.quadratures import GaussQuadrature1D

class TestShapeFunction(unittest.TestCase):
    def test_shape_function_value(self):
        shape_function = ShapeFunction1D()
        npt.assert_almost_equal(shape_function.value(0.0, 0), 0.5, 13)
        npt.assert_almost_equal(shape_function.value(0.0, 1), 0.5, 13)
        npt.assert_almost_equal(shape_function.value(1.0, 0), 0.0, 13)
        npt.assert_almost_equal(shape_function.value(-1.0, 0), 1.0, 13)
        npt.assert_almost_equal(shape_function.value(1.0, 1), 1.0, 13)
        npt.assert_almost_equal(shape_function.value(-1.0, 1), 0.0, 13)
        
    def test_shape_function_derivative(self):
        shape_function = ShapeFunction1D()
        npt.assert_almost_equal(shape_function.derivative(0.0, 0), -0.5, 13)
        npt.assert_almost_equal(shape_function.derivative(0.0, 1), 0.5, 13)
        npt.assert_almost_equal(shape_function.derivative(1.0, 0), -0.5, 13)
        npt.assert_almost_equal(shape_function.derivative(-1.0, 0), -0.5, 13)
        npt.assert_almost_equal(shape_function.derivative(1.0, 1), 0.5, 13)
        npt.assert_almost_equal(shape_function.derivative(-1.0, 1), 0.5, 13)
        
    def test_shape_function_integral(self):
        shape_function = ShapeFunction1D()
        quadrature = GaussQuadrature1D()
        func = lambda chi : shape_function.value(chi, 0)
        npt.assert_almost_equal(quadrature.integral(func), 1.0, 13)

        func = lambda chi : shape_function.value(chi, 1)
        npt.assert_almost_equal(quadrature.integral(func), 1.0, 13)

        func = lambda chi : shape_function.derivative(chi, 0)
        npt.assert_almost_equal(quadrature.integral(func), -1.0, 13)

        func = lambda chi : shape_function.derivative(chi, 1)
        npt.assert_almost_equal(quadrature.integral(func), 1.0, 13)


if __name__ == "__main__":
    unittest.main()
