__author__ = 'Vlad'

import unittest
import numpy as np
from main import dist_to_ellipse, mass_center, circle_fitting, ellipse_fitting, generate_points
from numpy import abs, sqrt
from numpy.linalg import norm

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

    def test_mass_center(self):
        prec = 1e-12
        center = [2, 3]
        radius = 4
        [C, R] = mass_center(generate_points(points_num=10, Xc=center, a=radius, b=radius, alpha=0))
        self.assertTrue(norm(C - center) < prec)
        self.assertTrue(abs(R - radius) < prec)

    def test_circle_fitting(self):
        # TODO: check the precision limit in a bigger test
        rel_prec = 1e-6
        center = [2, 3]
        radius = 4
        params = circle_fitting(generate_points(points_num=10, Xc=center, a=radius, b=radius, alpha=0), rel_prec)
        C = np.array([params[0], params[1]])
        R = params[2]
        self.assertTrue(norm(C - center)/norm(center) < rel_prec)
        self.assertTrue(abs(R - radius)/R < rel_prec)



if __name__ == '__main__':
    unittest.main()