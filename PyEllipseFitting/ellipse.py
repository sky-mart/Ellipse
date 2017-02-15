# _*_ coding: utf-8

__author__ = 'Vlad'

import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import solve, norm
from numpy import sin, cos, pi, arctan


# solve equation ax^2 + bx + c = 0
def solve2(a, b, c):
    if a == 0:
        if b == 0:
            return -1  # any x is a solution
        return [-c / b, -c / b]

    D = b*b - 4*a*c

    if D > 0:
        x1 = (-b + np.sqrt(D)) / (2*a)
        x2 = (-b - np.sqrt(D)) / (2*a)
        return [x1, x2]
    elif D < 0:
        x1 = complex(-b, np.sqrt(-D)) / (2*a)
        x2 = complex(-b, -np.sqrt(-D)) / (2*a)
        return [x1, x2]
    else:
        return [-b / (2*a), -b / (2*a)]


def mypow(x, p):
    if x >= 0:
        return x ** p
    return -((-x) ** p)


# solve equation ax^3 + bx^2 + cx + d = 0
def solve3(a, b, c, d):
    a = float(a); b = float(b); c = float(c); d = float(d)
    # reduce to the form y^3 + py + q = 0
    # substitution: x = y - b/3a
    p = c/a - b*b/(3*a*a)
    q = 2*b*b*b/(27*a*a*a) - b*c/(3*a*a) + d/a

    Q = p*p*p/27 + q*q/4
    # Q > 0 - one real root and two complex conjugated roots
    # Q = 0 - one single real root and one double real root, or,
	#         if p = q = 0, then one triple real root
	# Q < 0 - three real roots

    if Q >= 0:
        alpha   = mypow(-q/2 + np.sqrt(Q), 1.0/3)
        beta    = mypow(-q/2 - np.sqrt(Q), 1.0/3)
    else:
        alpha   = complex(-q/2, np.sqrt(-Q)) ** (1.0/3)
        beta    = complex(-q/2, -np.sqrt(-Q)) ** (1.0/3)

    x1 = alpha + beta - b/(3*a)
    x2 = complex(-(alpha+beta)/2, (alpha-beta)*np.sqrt(3)/2) - b/(3*a)
    x3 = complex(-(alpha+beta)/2, -(alpha-beta)*np.sqrt(3)/2) - b/(3*a)
    return [x1, x2, x3]


# solve equation ax^4 + bx^3 + cx^2 + dx + e = 0
def solve4(a, b, c, d, e):
    a = float(a); b = float(b); c = float(c); d = float(d); e = float(e)
    b /= a; c /= a; d /= a; e /= a;
    a = b; b = c; c = d; d = e;

    # reduce to the form y^4 + p*y^2 + q*y + r = 0
    p = b - 3*a*a/8
    q = a*a*a/8 - a*b/2 + c
    r = - 3*a*a*a*a/256 + a*a*b/16 - c*a/4 + d

    # obtain cubic resolvent A*s^3 + B*s^2 + C*s + D = 0
    A = 2
    B = -p
    C = -2*r
    D = r*p - q*q/4
    s1, s2, s3 = solve3(A, B, C, D)

    s = 0
    if np.real(s1) > p/2:
        s = s1
    elif np.real(s2) > p/2:
        s = s2
    elif np.real(s3) > p/2:
        s = s3

    a1 = 1; b1 = -np.sqrt(2*s-p); c1 = q/(2*np.sqrt(2*s-p)) + s
    a2 = 1; b2 = np.sqrt(2*s-p); c2 = -q/(2*np.sqrt(2*s-p)) + s

    x1, x2 = solve2(a1, b1, c1)
    x1 -= a/4
    x2 -= a/4
    x3, x4 = solve2(a2, b2, c2)
    x3 -= a/4
    x4 -= a/4
    return [x1, x2, x3, x4]


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
    roots = solve4(p[0], p[1], p[2], p[3], p[4])
    # roots = np.roots(p)
    lambdas = np.array([np.real(l) for l in roots if np.isreal(l)])
    distances = [np.sqrt((x[0] * l / (a ** 2 + l)) ** 2 + (x[1] * l / (b ** 2 + l)) ** 2) \
                 for l in lambdas if l != -a**2 and l != -b**2]
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


def ellipse_fitting_init_guess(points):
    Xc, _ = mass_center(points)

    closest = points[0]
    furthest = points[0]
    min_dist = norm(closest - Xc)
    max_dist = norm(furthest - Xc)

    for point in points[1:]:
        dist = norm(point - Xc)
        if dist < min_dist:
            min_dist = dist
            closest = point
        elif dist > max_dist:
            max_dist = dist
            furthest = point

    if abs(furthest[0] - Xc[0]) < 1e-6:
        alpha = np.pi / 2
    else:
        alpha = arctan((furthest[1] - Xc[1])/(furthest[0] - Xc[0]))
        # k1 = (furthest[1] - Xc[1])/(furthest[0] - Xc[0])
        # k2 = (-furthest[1]*Xc[0] + furthest[0]*Xc[1])/(furthest[0] - Xc[0])
    a = max_dist
    b = min_dist

    # initial conditions
    # TODO: check if there is singularity in case of circle
    # circle_Xc, circle_Yc, R = circle_fitting(global_points, rel_prec)
    # Xc = np.array([circle_Xc, circle_Yc])
    # a = R
    # b = R
    # alpha = 0

    # Xc, a, b, alpha = ellipse_fitting_init_guess(points)

    # k1, k2 = straight_line_lsfit(points)
    # np_k = np.polyfit(global_points[:, 0], global_points[:, 1], 1)
    # alpha = arctan(k1)
    # if dist_to_straight_line(k1, k2, closest) < dist_to_straight_line(k1, k2, furthest):
    #     a = min_dist
    #     b = max_dist
    # else:
    #     a = max_dist
    #     b = min_dist
    #
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
    # plt.plot(points[:, 0], points[:, 1], 'bo', \
    #          straight_x, straight_y, 'r', \
    #          small_circle_x, small_circle_y, 'g', \
    #          big_circle_x, big_circle_y, 'k')
    # plt.show()
    # visualize_fit(global_points, [Xc[0], Xc[1], a, b, alpha])

    return Xc, a, b, alpha


def ellipse_fitting(points, init, rel_prec=1e-6):
    N = len(points)
    J = np.zeros(shape=(N, 5), dtype=np.float32)  # jacobian

    Xc, a, b, alpha = init

    MAX_ITER_COUNT = 50
    MAX_ERROR_INCREASE_COUNT = 4
    iter_count = 1
    while iter_count <= MAX_ITER_COUNT:
        canon_points = np.array([global_to_canonical(X, alpha, Xc) for X in points])
        dls = [dist_to_ellipse(a, b, x) for x in canon_points]
        e = np.array([-item[0] for item in dls], dtype=np.float32)  # penalty


        if iter_count == 1:
            error_increase_count = 0
            prev_error = norm(e)
        else:
            prev_error = error
        error = norm(e)

        if error > prev_error:
            error_increase_count += 1
            if error_increase_count >= MAX_ERROR_INCREASE_COUNT:
                break
        else:
            error_increase_count = 0

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
        if abs(1 - abs(a/b)) < rel_prec:
            Js = np.dot(J[:, :4].transpose(), J[:, :4])
            es = np.dot(J[:, :4].transpose(), e)
            # da = solve(Js, es)
            da = linsolve(Js, es)
        else:
            Js = np.dot(J.transpose(), J)
            es = np.dot(J.transpose(), e)
            # da = solve(Js, es)
            da = linsolve(Js, es)
            alpha += da[4]

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

        if rel_prec_achieved or norm(da) < rel_prec**1.5:
            print iter_count
            return [Xc[0], Xc[1], abs(a), abs(b), alpha]

        iter_count += 1
    return [-1, -1, -1, -1, -1]


def generate_points(points_num, Xc, a, b, alpha, noise_level=0):
    step = 2 * np.pi / points_num
    points = np.zeros(shape=(points_num, 2), dtype=np.float32)

    t = 0
    for i in xrange(points_num):
        x = np.array([
            a * np.cos(t) * (1 + np.random.uniform(-noise_level, noise_level)),
            b * np.sin(t) * (1 + np.random.uniform(-noise_level, noise_level))
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


def ludcmp(A):
    N = A.shape[0]
    vv = np.zeros((N, 1))
    indx = np.zeros((N, 1), dtype=int)

    # calculate the biggest elements in every row
    for i in xrange(N):
        big = 0.0
        for j in xrange(N):
            temp = abs(A[i][j])
            if temp > big:
                big = temp
        if big == 0.0:
            print 'Singular matrix in routine ludcmp'
            return
        vv[i] = 1.0 / big

    for j in xrange(N):
        # calculate elements of L
        for i in xrange(j):
            sum = A[i, j]
            for k in xrange(i):
                sum -= A[i, k] * A[k, j]
            A[i, j] = sum

        # calculate elements of U
        big = 0.0
        for i in xrange(j, N):
            sum = A[i, j]
            for k in xrange(j):
                sum -= A[i, k] * A[k, j]
            A[i, j] = sum

            # find the row with the biggest jth element
            dum = vv[i] * abs(sum)
            if dum >= big:
                big = dum
                imax = i

        # figure out if we need to interchange rows
        if j != imax:
            for k in xrange(N):
                dum = A[imax, k]
                A[imax, k] = A[j, k]
                A[j, k] = dum
            vv[imax] = vv[j]
        # remember where the jth row goes
        indx[j] = imax

        # carry out division by the revealed biggest jj-th element
        dum = 1.0 / A[j, j]
        for i in xrange(j+1, N):
            A[i][j] *= dum

    return indx


def lubksb(A, b, indx):
    N = A.shape[0]
    ii = -1
    for i in xrange(N):
        ip = int(indx[i])
        sum = b[ip]
        b[ip] = b[i]
        if ii >= 0:
            for j in xrange(ii, i):
                sum -= A[i, j] * b[j]
        elif sum:
            ii = i
        b[i] = sum

    for i in np.arange(N-1, -1, -1):
        sum = b[i]
        for j in xrange(i+1, N):
            sum -= A[i, j] * b[j]
        b[i] = sum / A[i, i]
    return b

def linsolve(A, b):
    indx = ludcmp(A)
    return lubksb(A, b, indx)


if __name__ == '__main__':
    l = np.array([[1.0, 0.0, 0.0],
                  [1.0/7.0, 1.0, 0.0],
                  [0.0, 1.5, 1.0]])
    u = np.array([[7.0, 0.0, 2.0],
                  [0, 2.0, 33.0/7.0],
                  [0, 0, 13.5]])

    A = np.array([[1.0, 2.0, 5.0],
                  [0.0, 3.0, 21.0],
                  [7.0, 0.0, 2.0]])
    # print ludcmp(A)
    A_src = np.copy(A)
    #
    # N = A.shape[0]
    # LU, indx = ludcmp(A)

    # for i in xrange(N):
    #     for j in xrange(N):
    #         temp = A_src[i, j]
    #         A_src[i, j] = A_src[indx[i], j]
    #         A_src[indx[i], j] = temp
    #
    #
    # L = np.zeros((N, N))
    # U = np.zeros((N, N))
    # for i in xrange(N):
    #     L[i, i] = 1.0
    #     for j in xrange(i):
    #         L[i, j] = A[i, j]
    #     for j in xrange(i, N):
    #         U[i, j] = A[i, j]
    #
    # print LU
    # print L
    # print U
    # print np.dot(L, U)
    # print np.dot(L, U) - A_src


    b = np.array([3.0, 1.0, 2.5])
    b_src = np.copy(b)

    # for i in xrange(N):
    #     temp = b_src[i]
    #     b_src[i] = b_src[indx[i]]
    #     b_src[indx[i]] = temp

    # lubksb(A, b, indx)
    # print b
    print np.dot(A_src, linsolve(A, b)) - b_src

    # points = generate_points(500, [1.0, 1.0], 1.5, 1.0, 3.0543261909900767, noise_level=0.2)
    # init = ellipse_fitting_init_guess(points)
    # print ellipse_fitting(points, init)


