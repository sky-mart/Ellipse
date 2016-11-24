__author__ = 'Vlad'

import unittest
import numpy as np
from main import dist_to_ellipse
from numpy import abs, sqrt

class TestEllipse(unittest.TestCase):

    def test_dist_to_ellipse(self):
        prec = 1e-12
        self.assertTrue(abs(dist_to_ellipse(a=2, b=1, x=[4, 0])[0] - 2) < prec)
        self.assertTrue(abs(dist_to_ellipse(a=2, b=1, x=[0, 5])[0] - 4) < prec)
        self.assertTrue(abs(dist_to_ellipse(a=2, b=1, x=[-3, 0])[0] - 1) < prec)
        self.assertTrue(abs(dist_to_ellipse(a=2, b=1, x=[0, -8])[0] - 7) < prec)

        # circle tests
        self.assertTrue(abs(dist_to_ellipse(a=1, b=1, x=[2, 2])[0] - (sqrt(8) - 1)) < prec)
        self.assertTrue(abs(dist_to_ellipse(a=1, b=1, x=[-3, 4])[0] - (5 - 1)) < prec)

if __name__ == '__main__':
    unittest.main()