__author__ = 'Vlad'

import unittest
import numpy as np
from main import dist_to_ellipse, mass_center, circle_fitting, ellipse_fitting, generate_points, visualize_fit
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

    def angle_to_0_2pi(self, x, prec):
        while x < 0: #and abs(x) > prec:
            x += 2 * np.pi
        while x >= 2 * np.pi: #and abs(x - 2 * np.pi) > prec:
            x -= 2 * np.pi
        return x

    def compare_angles(self, alpha, beta, prec):
        if abs(alpha) < prec:
            alpha = 0
        if abs(beta) < prec:
            beta = 0
        alpha = self.angle_to_0_2pi(alpha, prec)
        beta = self.angle_to_0_2pi(beta, prec)
        return abs(alpha - beta) < prec

    def compare_ellipse_params(self, src_params, params, rel_prec):
        for i in xrange(2):
            if src_params[i] == 0:
                self.assertTrue(abs(params[i]) < rel_prec,
                                msg="src_params: " + str(src_params))
            else:
                self.assertTrue(abs(1 - params[i]/src_params[i]) < rel_prec,
                                msg="src_params: " + str(src_params))

        if abs(1 - params[2]/src_params[2]) < rel_prec and abs(1 - params[3]/src_params[3]) < rel_prec:
            if abs(1 - src_params[2]/src_params[3]) > rel_prec:
                self.assertTrue(self.compare_angles(params[4], src_params[4], rel_prec) or
                                self.compare_angles(params[4], src_params[4] + np.pi, rel_prec),
                                msg="src_params: " + str(src_params))
        elif abs(1 - params[3]/src_params[2]) < rel_prec and abs(1 - params[2]/src_params[3]) < rel_prec:
            if abs(1 - src_params[2]/src_params[3]) > rel_prec:
                self.assertTrue(self.compare_angles(params[4], src_params[4] + np.pi/2, rel_prec) or
                                self.compare_angles(params[4], src_params[4] + 3*np.pi/2, rel_prec),
                                msg="src_params: " + str(src_params))
        else:
            self.assertTrue(False, msg="src_params: " + str(src_params))

    def test_ellipse_fitting(self):
        rel_prec = 1e-5

        b = 1.0
        Xc = 1.0; Yc = 10;
        for points_num in xrange(30, 50):
            # for Xc in np.linspace(-1.0, 1.0, 10):
            #     for Yc in np.linspace(-1.0, 1.0, 10):
            for a in np.linspace(1.0, 2.0, 10):
                for alpha in np.linspace(0, 2*np.pi, 5):
                    src_params = [Xc, Yc, a, b, alpha]
                    points = generate_points(points_num, [Xc, Yc], a, b, alpha)

                    # visualize_fit(points, src_params)
                    params = ellipse_fitting(points, rel_prec)

                    self.compare_ellipse_params(src_params, params, rel_prec)
                    print points_num, a, alpha

if __name__ == '__main__':
    unittest.main()