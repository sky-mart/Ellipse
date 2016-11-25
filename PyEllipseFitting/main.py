__author__ = 'Vlad'

import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import solve, norm
from numpy import sin, cos, pi, fabs


def dist_to_ellipse(a, b, x):
    if x[0] == 0 and x[1] == 0:
        return min(a, b)

    # Equation from necessary condition of conditional extrema
    # Looking for point on ellipse (Xe, Ye):
    # (x-Xe)^2 + (y-Ye)^2 -> min, when Xe^2/a^2 + Ye^2/b^2 = 1
    p = [
        1,
        2 * (a ** 2 + b ** 2),
        a ** 4 + b ** 4 + 4 * a ** 2 * b ** 2 - a ** 2 * x[0] ** 2 - b ** 2 * x[1] ** 2,
        2.0 * a ** 2 * b ** 2 * (a ** 2 + b ** 2 - x[0] ** 2 - x[1] ** 2),
        a ** 2 * b ** 2 * (a ** 2 * b ** 2 - b ** 2 * x[0] ** 2 - a ** 2 * x[1] ** 2)
    ]

    # Choose lambda that gives minimal distance
    lambdas = np.array([np.real(l) for l in np.roots(p) if np.isreal(l)])
    distances = [np.sqrt((x[0] * l / (a ** 2 + l)) ** 2 + (x[1] * l / (b ** 2 + l)) ** 2) for l in lambdas]
    d = min(distances)
    l = lambdas[distances.index(d)]
    return d, l


def global_to_canonical(X, alpha, Xc):
    R = np.array([[cos(alpha), sin(alpha)],
                  [-sin(alpha), cos(alpha)]])

    return np.dot(R, X - Xc)


def canonical_to_global(x, alpha, Xc):
    R = np.array([[cos(alpha), -sin(alpha)],
                  [sin(alpha), cos(alpha)]])

    return np.dot(R, x) + Xc


def ellipse_fitting(global_points, rel_prec=1e-6):
    N = len(global_points)
    J = np.zeros(shape=(N, 5), dtype=np.float32)  # jacobian

    # initial conditions
    # TODO: check if there is singularity in case of circle
    circle_Xc, circle_Yc, R = circle_fitting(global_points, rel_prec)
    Xc = np.array([circle_Xc, circle_Yc])
    a = R
    b = R
    alpha = 0

    while True:
        canon_points = np.array([global_to_canonical(X, alpha, Xc) for X in global_points])
        dls = [dist_to_ellipse(a, b, x) for x in canon_points]
        e = np.array([-item[0] for item in dls], dtype=np.float32)  # penalty

        # fill jacobian
        for i in xrange(N):
            d, l = dls[i]
            x, y = canon_points[i]

            if l == 0 and d == 0:
                J[i] = [0, 0, 0, 0, 0]
                continue

            # auxiliary variables
            dg_dl = 4*l**3 + 6*l**2 * (a**2 + b**2) + \
                2*l * (a**4 + b**4 + 4*a**2*b**2 - a**2*x**2 - b**2*y**2) + \
                2*a**2*b**2 * (a**2 + b**2 - x**2 - y**2)

            dg_dx = -2*x*a**2 * (b**2 + l)**2

            dg_dy = -2*y*b**2 * (a**2 + l)**2

            dg_da = l**3 * 4*a + l**2 * 2*a * (2*a**2 + 4*b**2 - x**2) + \
                l * 4*a*b**2 * (2*a**2 + b**2 - x**2 - y**2) + \
                2*a*b**2 * (2*a**2*b**2 - b**2*x**2 - 2*a**2*y**2)

            dg_db = l**3 * 4*b + l**2 * 2*b * (2*b**2 + 4*a**2 - y**2) + \
                l * 4*a**2*b * (a**2 + 2*b**2 - x**2 - y**2) + \
                2*a**2*b * (2*a**2*b**2 - 2*b**2*x**2 - a**2*y**2)

            dl_dx = - dg_dx / dg_dl
            dl_dy = - dg_dy / dg_dl
            dl_da = - dg_da / dg_dl
            dl_db = - dg_db / dg_dl

            dx_dXc = -cos(alpha)
            dx_dYc = -sin(alpha)
            dx_dalpha = y
            dx_da = 0
            dx_db = 0

            dy_dXc = sin(alpha)
            dy_dYc = -cos(alpha)
            dy_dalpha = -x
            dy_da = 0
            dy_db = 0

            dl_dXc = dl_dx * dx_dXc + dl_dy * dy_dXc
            dl_dYc = dl_dx * dx_dYc + dl_dy * dy_dYc
            dl_dalpha = dl_dx * dx_dalpha + dl_dy * dy_dalpha

            da2_dXc = 0
            da2_dYc = 0
            da2_da = 2*a
            da2_db = 0
            da2_dalpha = 0

            db2_dXc = 0
            db2_dYc = 0
            db2_da = 0
            db2_db = 2*b
            db2_dalpha = 0

            dd_dXc = (l/d) * (
                x / (a**2 + l)**3 * (a**2*x * dl_dXc + l*(a**2 + l) * dx_dXc - l*x * da2_dXc)
                +
                y / (b**2 + l)**3 * (b**2*y * dl_dXc + l*(b**2 + l) * dy_dXc - l*y * db2_dXc)
            )

            dd_dYc = (l/d) * (
                x / (a**2 + l)**3 * (a**2*x * dl_dYc + l*(a**2 + l) * dx_dYc - l*x * da2_dYc)
                +
                y / (b**2 + l)**3 * (b**2*y * dl_dYc + l*(b**2 + l) * dy_dYc - l*y * db2_dYc)
            )

            dd_da = (l/d) * (
                x / (a**2 + l)**3 * (a**2*x * dl_da + l*(a**2 + l) * dx_da - l*x * da2_da)
                +
                y / (b**2 + l)**3 * (b**2*y * dl_da + l*(b**2 + l) * dy_da - l*y * db2_da)
            )

            dd_db = (l/d) * (
                x / (a**2 + l)**3 * (a**2*x * dl_db + l*(a**2 + l) * dx_db - l*x * da2_db)
                +
                y / (b**2 + l)**3 * (b**2*y * dl_db + l*(b**2 + l) * dy_db - l*y * db2_db)
            )

            dd_dalpha = (l/d) * (
                x / (a**2 + l)**3 * (a**2*x * dl_dalpha + l*(a**2 + l) * dx_dalpha - l*x * da2_dalpha)
                +
                y / (b**2 + l)**3 * (b**2*y * dl_dalpha + l*(b**2 + l) * dy_dalpha - l*y * db2_dalpha)
            )

            J[i] = np.array([
                dd_dXc,
                dd_dYc,
                dd_da,
                dd_db,
                dd_dalpha
            ])

        # solve linear system J * da = e
        # to find ellipse parameters' changes

        if abs(1 - a/b) < rel_prec:
            Js = np.dot(J[:, :4].transpose(), J[:, :4])
            es = np.dot(J[:, :4].transpose(), e)
            da = solve(Js, es)
        else:
            Js = np.dot(J.transpose(), J)
            es = np.dot(J.transpose(), e)
            da = solve(Js, es)
            alpha += da[4]

        Xc[0] += da[0]
        Xc[1] += da[1]
        a += da[2]
        b += da[3]
        params = [Xc[0], Xc[1], a, b, alpha]

        print "norm:", norm(da)
        rel_prec_achieved = True
        for i in xrange(len(da)):
            if abs(da[i] / params[i]) > rel_prec:
                rel_prec_achieved = False
                break
        if rel_prec_achieved or norm(da) < 1e-12:
            return [Xc[0], Xc[1], a, b, alpha]


def generate_points(points_num, Xc, a, b, alpha, noise_level=0):
    step = 2 * np.pi / points_num
    points = np.zeros(shape=(points_num, 2), dtype=np.float32)

    t = 0
    for i in xrange(points_num):
        x = np.array([
            a * np.cos(t),
            b * np.sin(t)
        ])
        points[i] = canonical_to_global(x, -alpha, Xc)
        t += step

    return points


def circle_fitting(points, rel_prec=1e-6):
    N = len(points)
    J = np.zeros(shape=(N, 3), dtype=np.float32)
    e = np.zeros(shape=(N, 1), dtype=np.float32)

    Xc, R = mass_center(points)

    while True:
        for i, X in enumerate(points):
            Dc = norm(X - Xc)
            e[i] = -abs(Dc - R)

            if Dc == 0:
                J[i, 0] = 0
                J[i, 1] = 0
                J[i, 2] = 1
            elif Dc > R:
                J[i, 0] = (Xc[0] - X[0]) / Dc
                J[i, 1] = (Xc[1] - X[1]) / Dc
                J[i, 2] = -1
            else:
                J[i, 0] = (X[0] - Xc[0]) / Dc
                J[i, 1] = (X[1] - Xc[1]) / Dc
                J[i, 2] = 1

        Js = np.dot(J.transpose(), J)
        es = np.dot(J.transpose(), e)
        da = solve(Js, es)

        Xc[0] += da[0]
        Xc[1] += da[1]
        R += da[2]

        # print "norm:", norm(da)
        if abs(da[0])/Xc[0] < rel_prec and abs(da[1])/Xc[1] < rel_prec and abs(da[2])/R:
            return [np.asscalar(Xc[0]), np.asscalar(Xc[1]), np.asscalar(R)]


def mass_center(points):
    C = sum(points) / len(points)
    R = np.array([norm(x - C) for x in points]).mean()
    return C, R

points = generate_points(30, [-1, -1], 2.2, 1.5, 2*pi/9)
params = ellipse_fitting(points)
print params
plt.plot(points[:, 0], points[:, 1], 'bo')
plt.show()


# dd_dXc = (l/d) * (
#     a**2*x**2 / (a**2 + l)**3 *
#     (cos(alpha) * (dg_dx/dg_dl + l*(a**2 + l)/a**2/x) - sin(alpha) * dg_dy/dg_dl)
#     +
#     b**2*y**2 / (b**2 + l)**3 *
#     (-sin(alpha) * (dg_dy/dg_dl + l*(b**2 + l)/b**2/y) + cos(alpha) * dg_dx/dg_dl)
# )
#
# dd_dYc = (l / d) * (
#     a**2*x**2 / (a**2 + l)**3 *
#     (sin(alpha) * (dg_dx/dg_dl + l*(a**2 + l)/a**2/x) + cos(alpha) * dg_dy/dg_dl)
#     +
#     b**2*y**2 / (b**2 + l)**3 *
#     (cos(alpha) * (dg_dy/dg_dl + l*(b**2 + l)/b**2/y) + sin(alpha) * dg_dx/dg_dl)
# )
#
# dd_da = -(l/d) * (
#     a*x**2 * (a*dg_da/dg_dl + 2*l) / (a**2 + l)**3
#     +
#     b**2*x**2 * dg_da/dg_dl / (b**2 + l)**3
# )
#
# dd_db = -(l/d) * (
#     a**2*x**2 * dg_db/dg_dl / (a**2 + l)**3
#     +
#     b*x**2 * (b*dg_db/dg_dl + 2*l) / (b**2 + l)**3
# )
#
# dd_dalpha = (l/d) * (
#     a**2*x**2 / (a**2 + l)**3 *
#     (-y * (dg_dx/dg_dl + l*(a**2 + l)/a**2/x) - x * dg_dy/dg_dl)
#     +
#     b**2*y**2 / (b ** 2 + l) ** 3 *
#     (x * (dg_dy/dg_dl + l*(b**2 + l)/b**2/y) - y * dg_dx/dg_dl)
# )

