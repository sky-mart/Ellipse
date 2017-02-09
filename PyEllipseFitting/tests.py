__author__ = 'Vlad'

from nose_parameterized import parameterized

import unittest
import numpy as np
import ellipse
from ellipse import dist_to_ellipse, mass_center, circle_fitting, \
    ellipse_fitting, ellipse_fitting_init_guess, generate_points, visualize_fit
from numpy import sqrt
from numpy.linalg import norm


def cartesian(arrays, out=None):
    """
    Generate a cartesian product of input arrays.

    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the cartesian product of.
    out : ndarray
        Array to place the cartesian product in.

    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.

    Examples
    --------
    >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],
           [3, 5, 6],
           [3, 5, 7]])

    """

    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = n / arrays[0].size
    out[:, 0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m, 1:])
        for j in xrange(1, arrays[0].size):
            out[j * m:(j + 1) * m, 1:] = out[0:m, 1:]
    return out


class TestEllipse(unittest.TestCase):
    @parameterized.expand([
        (1, 5, 6, [-2, -3]),                    # D > 0
        (1, 1, 2.5, [-0.5+1.5j, -0.5-1.5j]),    # D < 0
        (1, -4, 4, [2, 2]),                     # D == 0
        (0, 5, 5, [-1, -1]),                    # a == 0, b != 0, c != 0
        (0, 0, 3, -1),                          # a == 0, b == 0, c != 0
        (0, 3, 0, [0, 0]),                      # a == 0, b != 0, c == 0
        (1, 5, 0, [0, -5]),                     # a != 0, b != 0, c == 0
        (1, 0, 9, [3j, -3j])                    # a != 0, b == 0, c != 0
    ])
    def test_solve2(self, a, b, c, roots):
        self.assertTrue(ellipse.solve2(a, b, c) == roots)

    @parameterized.expand([
        (1, 5, 6, -3),                          # Q > 0
        (1, 3, -3, -1),                         # Q < 0
        (1, 6, 9, -1),                          # Q == 0
        (1, 6, 9, 22),                          # random
    ])
    def test_solve3(self, a, b, c, d):
        prec = 1e-6
        roots = ellipse.solve3(a, b, c, d)
        for r in roots:
            eps = abs(a*r*r*r + b*r*r + c*r + d)
            self.assertTrue(eps < prec)

    @parameterized.expand([
        (1, 5, 6, -3, -21)
    ])
    def test_solve4(self, a, b, c, d, e):
        prec = 1e-6
        roots = ellipse.solve4(a, b, c, d, e)
        for r in roots:
            eps = abs(a*r*r*r*r + b*r*r*r + c*r*r + d*r + e)
            self.assertTrue(eps < prec)

    @parameterized.expand([
        (2, 1, [4, 0], 2),
        (2, 1, [0, 5], 4),
        (2, 1, [-3, 0], 1),
        (2, 1, [0, -8], 7),
        (1, 1, [2, 2], sqrt(8)-1),
        (1, 1, [-3, 4], 5-1)
    ])
    def test_dist_to_ellipse(self, a, b, x, d):
        prec = 1e-12
        self.assertTrue(abs(dist_to_ellipse(a, b, x)[0] - d) < prec)
        # self.assertTrue(abs(dist_to_ellipse(a=2, b=1, x=[0, 5])[0] - 4) < prec)
        # self.assertTrue(abs(dist_to_ellipse(a=2, b=1, x=[-3, 0])[0] - 1) < prec)
        # self.assertTrue(abs(dist_to_ellipse(a=2, b=1, x=[0, -8])[0] - 7) < prec)
        #
        # # circle tests
        # self.assertTrue(abs(dist_to_ellipse(a=1, b=1, x=[2, 2])[0] - (sqrt(8) - 1)) < prec)
        # self.assertTrue(abs(dist_to_ellipse(a=1, b=1, x=[-3, 4])[0] - (5 - 1)) < prec)

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
        self.assertTrue(norm(C - center) / norm(center) < rel_prec)
        self.assertTrue(abs(R - radius) / R < rel_prec)

    # @parameterized.expand(cartesian((
    #     xrange(30, 50),
    #     np.linspace(1.0, 2.0, 10),
    #     np.linspace(0, 2*np.pi, 5)))
    # )
    # def test_ellipse_fitting(self, points_num, a, alpha):
    #     rel_prec = 1e-6
    #     Xc = 1.0
    #     Yc = 10
    #     b = 1.0
    #
    #     src_params = [Xc, Yc, a, b, alpha]
    #     points = generate_points(points_num, [Xc, Yc], a, b, alpha)
    #     params = ellipse_fitting(points, rel_prec)
    #     self.compare_ellipse_params(src_params, params, rel_prec)

    @parameterized.expand(cartesian((
        np.linspace(500, 501, 1),
        np.linspace(1.0, 3.0, 4),
        np.linspace(-np.pi+np.pi/36, np.pi-np.pi/36, 4),
        np.linspace(0.2, 0.3, 4)))
    )
    def test_ellipse_fitting_noisy(self, points_num, a, alpha, noise_level):
        rel_prec = 3e-1
        Xc = 1.0
        Yc = 1.0
        b = 1.0

        src_params = [Xc, Yc, a, b, alpha]
        points = generate_points(int(points_num), [Xc, Yc], a, b, alpha, noise_level)
        init = ellipse_fitting_init_guess(points)
        params = ellipse_fitting(points, init)
        if params == [-1, -1, -1, -1, -1]:
            self.assertTrue(False, msg='rashod: ' + "src_params: " + str(src_params))
        else:
            self.compare_ellipse_params(src_params, params, rel_prec)

    def compare_ellipse_params(self, src_params, params, rel_prec):
        msg = "src_params: " + str(src_params) + "\nparams: " + str(params)
        for i in xrange(2):
            if src_params[i] == 0:
                self.assertTrue(abs(params[i]) < rel_prec, msg=msg)
            else:
                self.assertTrue(abs(1 - params[i] / src_params[i]) < rel_prec, msg=msg)

        if abs(1 - params[2] / src_params[2]) < rel_prec and abs(1 - params[3] / src_params[3]) < rel_prec:
            if abs(1 - src_params[2] / src_params[3]) > rel_prec:
                self.assertTrue(self.compare_angles(params[4], src_params[4], rel_prec) or
                                self.compare_angles(params[4], src_params[4] + np.pi, rel_prec),
                                msg=msg)
        elif abs(1 - params[3] / src_params[2]) < rel_prec and abs(1 - params[2] / src_params[3]) < rel_prec:
            if abs(1 - src_params[2] / src_params[3]) > rel_prec:
                self.assertTrue(self.compare_angles(params[4], src_params[4] + np.pi / 2, rel_prec) or
                                self.compare_angles(params[4], src_params[4] + 3 * np.pi / 2, rel_prec),
                                msg=msg)
        else:
            self.assertTrue(False, msg=msg)

    def compare_angles(self, alpha, beta, prec):
        if abs(alpha) < 1e-6:
            alpha = 0
        if abs(beta) < 1e-6:
            beta = 0
        alpha = self.angle_to_0_2pi(alpha, prec)
        beta = self.angle_to_0_2pi(beta, prec)
        return abs(alpha - beta) < prec

    def angle_to_0_2pi(self, x, prec):
        while x < 0: # and abs(x) > prec:
            x += 2 * np.pi
        while x >= 2 * np.pi: # and abs(x - 2 * np.pi) > prec:
            x -= 2 * np.pi
        return x


if __name__ == '__main__':
    unittest.main()
    # print cartesian([[1, 2], [3, 4], [5, 6]])