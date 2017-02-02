__author__ = 'Vlad'

import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import solve, norm
from numpy import sin, cos, pi, arctan


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

def dist_to_straight_line(k1, k2, x):
    return abs(x[1] - k1*x[0] - k2) / k1

def ellipse_fitting(global_points, rel_prec=1e-6):
    N = len(global_points)
    J = np.zeros(shape=(N, 5), dtype=np.float32)  # jacobian

    # initial conditions
    # TODO: check if there is singularity in case of circle
    # circle_Xc, circle_Yc, R = circle_fitting(global_points, rel_prec)
    # Xc = np.array([circle_Xc, circle_Yc])
    # a = R
    # b = R
    # alpha = 0
    Xc, _ = mass_center(global_points)
    min_dist = norm(global_points[0] - Xc)
    closest = global_points[0]
    max_dist = norm(global_points[0] - Xc)
    furthest = global_points[0]
    for point in global_points[1:]:
        dist = norm(point - Xc)
        if dist < min_dist:
            min_dist = dist
            closest = point
        elif dist > max_dist:
            max_dist = dist
            furthest = point

    alpha = arctan((furthest[1] - Xc[1])/(furthest[0] - Xc[0]))
    a = max_dist
    b = min_dist
    # k1, k2 = straight_line_lsfit(global_points)
    # # np_k = np.polyfit(global_points[:, 0], global_points[:, 1], 1)
    # alpha = arctan(k1)
    # # TODO: check if x axis value is always bigger
    # if dist_to_straight_line(k1, k2, closest) < dist_to_straight_line(k1, k2, furthest):
    #     a = min_dist
    #     b = max_dist
    # else:
    #     a = max_dist
    #     b = min_dist

    # straight_x = np.linspace(Xc[0]-max_dist, Xc[0]+max_dist)
    # straight_y = straight_x * k1 + k2
    #
    # t = np.linspace(0, 2*np.pi, 50)
    # small_circle_x = Xc[0] + min_dist*cos(t)
    # small_circle_y = Xc[1] + min_dist*sin(t)
    #
    # big_circle_x = Xc[0] + max_dist*cos(t)
    # big_circle_y = Xc[1] + max_dist*sin(t)
    #
    # plt.plot(global_points[:, 0], global_points[:, 1], 'bo', \
    #          straight_x, straight_y, 'r', \
    #          small_circle_x, small_circle_y, 'g', \
    #          big_circle_x, big_circle_y, 'k')
    # plt.show()
    # visualize_fit(global_points, [Xc[0], Xc[1], a, b, alpha])

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
        # print "error:", norm(e) ** 2

        # solve linear system J * da = e
        # to find ellipse parameters' changes

        # try:
        if abs(1 - abs(a/b)) < rel_prec:
            Js = np.dot(J[:, :4].transpose(), J[:, :4])
            es = np.dot(J[:, :4].transpose(), e)
            da = solve(Js, es)
        else:
            Js = np.dot(J.transpose(), J)
            es = np.dot(J.transpose(), e)
            da = solve(Js, es)
            alpha += da[4]
        # except np.linalg.LinAlgError as ex:
        #     print ex

        Xc[0] += da[0]
        Xc[1] += da[1]
        a += da[2]
        b += da[3]
        params = [Xc[0], Xc[1], a, b, alpha]

        # visualize_fit(points, params)
        # print "error", norm(e)
        # print params
        # print da
        # print '\n'

        Xc_ok = [False, False]
        for i in xrange(2):
            Xc_ok[i] = check_pot_zero_param(Xc[i], da[i], rel_prec)

        alpha_ok = False
        if len(da) < 5:
            alpha_ok = True
        elif check_pot_zero_param(alpha, da[4], rel_prec):
            alpha_ok = True

        rel_prec_achieved = Xc_ok[0] and Xc_ok[1] and alpha_ok and \
                            abs(da[2]/a) < rel_prec and abs(da[3]/b) < rel_prec

        if rel_prec_achieved or norm(da) < rel_prec**2:
            return [Xc[0], Xc[1], abs(a), abs(b), alpha]


def generate_points(points_num, Xc, a, b, alpha, noise_level=0):
    step = 2 * np.pi / points_num
    points = np.zeros(shape=(points_num, 2), dtype=np.float32)

    t = 0
    for i in xrange(points_num):
        x = np.array([
            a * np.cos(t),
            b * np.sin(t)
        ])
        points[i] = canonical_to_global(x, alpha, Xc)
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
        Xc_ok = [False, False]
        for i in xrange(2):
            Xc_ok[i] = check_pot_zero_param(Xc[i], da[i], rel_prec)

        if Xc_ok[0] and Xc_ok[1] and abs(da[2])/R:
            return [np.asscalar(Xc[0]), np.asscalar(Xc[1]), np.asscalar(R)]


# check potentially zero parameter
def check_pot_zero_param(value, delta, rel_prec):
    if abs(value) < rel_prec:  # == 0
        return True if abs(delta) < rel_prec else False
    else:
        return True if abs(delta/value) < rel_prec else False


def mass_center(points):
    C = sum(points) / len(points)
    R = np.array([norm(x - C) for x in points]).mean()
    return C, R


def visualize_fit(points, params):
    fitted = generate_points(30, [params[0], params[1]], params[2], params[3], params[4])
    plt.plot(points[:, 0], points[:, 1], 'bo', fitted[:, 0], fitted[:, 1])
    plt.show()

def straight_line_lsfit(points):
    sx = 0; sy = 0; sx2 = 0; sxy = 0;
    for p in points:
        sx += p[0]
        sy += p[1]
        sx2 += p[0] * p[0]
        sxy += p[0] * p[1]

    a11 = sx;   a12 = len(points);
    a21 = sx2;  a22 = sx;
    b1 = sy;    b2 = sxy;

    k = (b1 * a22 - a12 * b2) / (a11 * a22 - a12 * a21)
    b = (a11 * b2 - b1 * a21) / (a11 * a22 - a12 * a21)
    return k, b


# points = generate_points(30, (-1.0, -1.0), 1.0, 1.6363636363636365, 10*np.pi/9)
# # # # plt.plot(points[:, 0], points[:, 1], 'bo')
# # # # plt.show()
# params = ellipse_fitting(points)
# print params
# visualize_fit(points, params)


# First error type: singular matrix
# 30 [-1.0, 0] 3 1.2 2*pi/9
# 30 [-1.0, 0] 3 1.3 2*pi/9

# Second error type: program does not respond
# 30 [-1.0, 0] 3 1.4 2*pi/9

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

