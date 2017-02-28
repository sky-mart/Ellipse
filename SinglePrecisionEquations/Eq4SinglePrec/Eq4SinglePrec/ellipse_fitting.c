//
//  ellipse_fitting.c
//  Eq4SinglePrec
//
//  Created by Влад Агиевич on 28.02.17.
//  Copyright © 2017 Mocsmart. All rights reserved.
//

#include "ellipse_fitting.h"

#include <math.h>
#include <float.h>      // FLT_MAX
#include <stdlib.h>     // malloc
#include <string.h>     // memset
#include <stdio.h>

void cparts(complexr z, real *pr, real *pi)
{
    *pr = crealr(z);
    *pi = crealr(z);
}

void cprint(complexr z)
{
    printf("%e + %ei", crealr(z), cimagr(z));
}

//real mypowr(real b, real e)
//{
//    return (b >= 0) ? powr(b, e) : -powr(-b, e);
//}

/*
 * solves equation ax^2 + bx + c = 0
 * puts roots in the array
 * returns 1 on success
 */
int solve2(real a, real b, real c, complexr roots[2])
{
    real D;
    if (a == 0) {
        if (b == 0) {
            return 0; // any x is a solution
        }
        roots[0] = roots[1] = -c / b;
        return 1;
    }
    
    D = b*b - 4*a*c;
    if (D >= 0) {
        roots[0] = (-b + sqrtr(D)) / (2*a);
        roots[1] = (-b - sqrtr(D)) / (2*a);
    } else {
        roots[0] = (-b + sqrtr(-D) * I) / (2*a);
        roots[1] = (-b - sqrtr(-D) * I) / (2*a);
    }
    return 1;
}

/*
 * solves equation ax^3 + bx^2 + cx + d = 0
 * puts roots in the array
 * returns 1 on success
 */
int solve3(real a, real b, real c, real d, complexr roots[3])
{
    real p, q, Q;
    complexr alpha, beta;
    //    real rp[5], ip[5]; // debug
#ifdef SINGLE_PREC_DEBUG
    printf("solve3:\n");
    printf("a=%e\n", a);
    printf("b=%e\n", b);
    printf("c=%e\n", c);
    printf("d=%e\n", d);
    printf("\n");
#endif
    // reduce to the form y^3 + py + q = 0
    // substitution: x = y - b/3a
    p = c/a - b*b/(3*a*a);
    q = 2*b*b*b/(27*a*a*a) - b*c/(3*a*a) + d/a;
    
    Q = p*p*p/27 + q*q/4;
#ifdef SINGLE_PREC_DEBUG
    printf("p=%e\n", p);
    printf("q=%e\n", q);
    printf("Q=%e\n", Q);
    printf("\n");
#endif
    // Q > 0 - one real root and two complex conjugated roots
    // Q = 0 - one single real root and one double real root, or,
    //         if p = q = 0, then one triple real root
    // Q < 0 - three real roots
    
    if (Q >= 0) {
        //alpha   = mypowr(-q/2 + sqrtr(Q), 1.0/3);
        //beta    = mypowr(-q/2 - sqrtr(Q), 1.0/3);
        alpha   = cbrtr(-q/2 + sqrtr(Q));
        beta    = cbrtr(-q/2 - sqrtr(Q));
    } else {
        alpha   = cpowr(-q/2 + sqrtr(-Q) * I, 1.0/3);
        beta    = cpowr(-q/2 - sqrtr(-Q) * I, 1.0/3);
    }
#ifdef SINGLE_PREC_DEBUG
    printf("alpha="); cprint(alpha); printf("\n");
    printf("beta="); cprint(beta); printf("\n");
    printf("\n");
#endif
    //    cparts(alpha, &rp[0], &ip[0]);
    //    cparts(beta, &rp[1], &ip[1]);
    
    roots[0] = alpha + beta - b/(3*a);
    roots[1] = -(alpha+beta)/2 - b/(3*a) + (alpha-beta)*sqrtr(3)/2 * I;
    roots[2] = -(alpha+beta)/2 - b/(3*a) - (alpha-beta)*sqrtr(3)/2 * I;
#ifdef SINGLE_PREC_DEBUG
    printf("roots[0]="); cprint(roots[0]); printf("\n");
    printf("roots[1]="); cprint(roots[1]); printf("\n");
    printf("roots[2]="); cprint(roots[2]); printf("\n");
    printf("\n\n");
#endif
    
    //    cparts(roots[0], &rp[2], &ip[2]);
    //    cparts(roots[1], &rp[3], &ip[3]);
    //    cparts(roots[2], &rp[4], &ip[4]);
    
    return 1;
}

/*
 * solves equation ax^4 + bx^3 + cx^2 + dx + e = 0
 * puts roots in the array
 * returns 1 on success
 */
int solve4(real a, real b, real c, real d, real e, complexr roots[4])
{
    real p, q, r;
    real A, B, C, D, s;
    complexr cube[3];
    real a1, b1, c1, a2, b2, c2;
    int i;
    //real rp[3], ip[3]; // debug
#ifdef SINGLE_PREC_DEBUG
    printf("solve4:\n");
    printf("a=%e\n", a);
    printf("b=%e\n", b);
    printf("c=%e\n", c);
    printf("d=%e\n", d);
    printf("e=%e\n", e);
    printf("\n");
#endif
    
    b /= a; c /= a; d /= a; e /= a;
    a = b; b = c; c = d; d = e;
    
    // reduce to the form y^4 + p*y^2 + q*y + r = 0
    p = b - 3*a*a/8;
    q = a*a*a/8 - a*b/2 + c;
    r = - 3*a*a*a*a/256 + a*a*b/16 - c*a/4 + d;
#ifdef SINGLE_PREC_DEBUG
    printf("p=%e\n", p);
    printf("q=%e\n", q);
    printf("r=%e\n", r);
    printf("\n");
#endif
    // obtain cubic resolvent A*s^3 + B*s^2 + C*s + D = 0
    A = 2;
    B = -p;
    C = -2*r;
    D = r*p - q*q/4;
    if (!solve3(A, B, C, D, cube)) {
        return 0;
    }
    
    //    for (i = 0; i < 3; i++) {
    //        rp[i] = crealr(cube[i]);
    //        ip[i] = cimagr(cube[i]);
    //    }
    
    s = 0;
    if (crealr(cube[0]) > p/2) {
        s = crealr(cube[0]);
    } else if (crealr(cube[1]) > p/2) {
        s = crealr(cube[1]);
    } else if (crealr(cube[2]) > p/2) {
        s = crealr(cube[2]);
    }
#ifdef SINGLE_PREC_DEBUG
    printf("s=%e\n", s);
    printf("\n");
#endif
    a1 = 1; b1 = -sqrtr(2*s-p); c1 = q/(2*sqrtr(2*s-p)) + s;
    a2 = 1; b2 = sqrtr(2*s-p); c2 = -q/(2*sqrtr(2*s-p)) + s;
#ifdef SINGLE_PREC_DEBUG
    printf("a1=%e\n", a1);
    printf("b1=%e\n", b1);
    printf("c1=%e\n", c1);
    printf("a2=%e\n", a2);
    printf("b2=%e\n", b2);
    printf("c2=%e\n", c2);
    printf("\n");
#endif
    
    if (!solve2(a1, b1, c1, &roots[0]) || !solve2(a2, b2, c2, &roots[2])) {
        return 0;
    }
    for (i = 0; i < 4; i++) {
        roots[i] -= a/4;
    }
#ifdef SINGLE_PREC_DEBUG
    printf("roots[0]="); cprint(roots[0]); printf("\n");
    printf("roots[1]="); cprint(roots[1]); printf("\n");
    printf("roots[2]="); cprint(roots[2]); printf("\n");
    printf("roots[3]="); cprint(roots[3]); printf("\n");
    printf("\n\n");
#endif
    return 1;
}

/*
 * Calculates the distance from the point x to the ellipse with axis a, b
 * *pd - output distance
 * *pl - output auxilliary parameter lambda
 * Returns 1 on success
 */
int dist_to_ellipse(real a, real b, real x[2], real *pd, real *pl)
{
    complexr roots[4];
    real d; // distance
    real l; // lambda
    int i;
    
    if (x[0] == 0 && x[1] == 0)
        return fminr(a, b);
#ifdef SINGLE_PREC_DEBUG
    printf("dist_to_ellipse:\n");
    printf("a=%e\n", a);
    printf("b=%e\n", b);
    printf("x0=%e\n", x[0]);
    printf("x1=%e\n", x[1]);
    printf("\n\n");
#endif
    // Equation from necessary condition of conditional extrema
    // Looking for point on ellipse (Xe, Ye):
    // (x-Xe)^2 + (y-Ye)^2 -> min, when Xe^2/a^2 + Ye^2/b^2 = 1
    if (!solve4(1.0f,
                2.0f * (a*a + b*b),
                a*a*a*a + b*b*b*b + 4*a*a*b*b - a*a*x[0]*x[0] - b*b*x[1]*x[1],
                2.0f * a*a*b*b * (a*a + b*b - x[0]*x[0] - x[1]*x[1]),
                a*a*b*b * (a*a*b*b - b*b*x[0]*x[0] - a*a*x[1]*x[1]),
                roots)) {
        return 0;
    }
    //    for (i = 0; i < 4; i++) {
    //        printf("r[%i] = ", i);
    //        cprint(roots[i]);
    //        printf("\n");
    //    }
    
    // Choose lambda that gives minimal distance
    *pd = FLT_MAX;
    for (i = 0; i < 4; i++) {
        if (cimagr(roots[i]) == 0) {
            l = crealr(roots[i]);
            
#define DIST_TO_ELLIPSE_EQ_PREC 1e-3f
            if ((fabsr(a*a + l) > DIST_TO_ELLIPSE_EQ_PREC) &&
                (fabsr(b*b + l) > DIST_TO_ELLIPSE_EQ_PREC))
            {
                d = sqrtr(powr(x[0]*l/(a*a + l), 2) + powr(x[1]*l/(b*b + l), 2));
                if (d < *pd) {
                    *pd = d;
                    *pl = l;
                }
            }
        }
    }
    return 1;
}

void dotprodr(real *a, real *b, uint n, real *pResult)
{
    uint i;
    *pResult = 0;
    for (i = 0; i < n; i++) {
        *pResult += a[i] * b[i];
    }
}

void vec_sub(real *a, real *b, real *d, uint n)
{
    uint i;
    for (i = 0; i < n; i++) {
        d[i] = a[i] - b[i];
    }
}

void vec_max(real *v, uint n, real *pResult, uint *pIndex)
{
    uint i, bigIndex = 0;
    real big = v[0];
    for (i = 1; i < n; i++) {
        if (v[i] > big) {
            bigIndex = i;
            big = v[i];
        }
    }
    *pResult = big;
    *pIndex = bigIndex;
}

void mat_init(matrix *pA, uint16_t nRows, uint16_t nCols, real *pData)
{
    pA->numRows = nRows;
    pA->numCols = nCols;
    pA->pData   = pData;
}

int mat_trans(matrix *pA, matrix *pAt)
{
    int i, j;
    for (i = 0; i < pA->numRows; i++) {
        for (j = 0; j < pA->numCols; j++) {
            pAt->pData[j * pAt->numCols + i] = pA->pData[i * pA->numCols + j];
        }
    }
    return 1;
}

int mat_sub(matrix *pA, matrix *pB, matrix *pC)
{
    int i, j;
    for (i = 0; i < pA->numRows; i++) {
        for (j = 0; j < pA->numCols; j++) {
            pC->pData[i * pC->numCols + j] = pA->pData[i * pA->numCols + j] -
                                                pB->pData[i *pB->numCols + j];
        }
    }
    return 1;
}

int mat_mult(matrix *pA, matrix *pB, matrix *pC)
{
    int i, j, k;
    real tmp;
    for (i = 0; i < pC->numRows; i++) {
        for (j = 0; j < pC->numCols; j++) {
            tmp = 0;
            for (k = 0; k < pA->numCols; k++) {
                tmp += pA->pData[i * pA->numCols + k] * pB->pData[k * pB->numCols + j];
            }
            pC->pData[i * pC->numCols + j] = tmp;
        }
    }
    return 1;
}

real norm(real *v, uint n)
{
    real result;
    dotprod(v, v, n, &result);
    return sqrtr(result);
}

real dist(real *u, real *v, int n)
{
    real val;
    real *d = (real *)malloc(n * sizeof(real));
    
    vsub(u, v, d, n);
    val = norm(d, n);
    free(d);
    return val;
}

/*
 * Initial guess for ellipse_fitting routine
 * Returns Xc, a, b, alpha in pointers
 */
void ellipse_fitting_init_guess(real **points, int N, struct ellipse *pInit)
{
    real R;
    int closest, furthest, i;
    real min_dist, max_dist, d;
    mass_center(points, N, pInit->Xc, &R);
    
    closest = furthest = 0;
    min_dist = dist(points[closest], pInit->Xc, 2);
    max_dist = dist(points[furthest], pInit->Xc, 2);
    
    for (i = 1; i < N; i++) {
        d = dist(points[i], pInit->Xc, 2);
        if (d < min_dist) {
            min_dist = d;
            closest = i;
        } else if (d > max_dist) {
            max_dist = d;
            furthest = i;
        }
    }
    
    if (fabsr(points[furthest][0] - pInit->Xc[0]) < 1e-6) {
        pInit->alpha = M_PI / 2;
    } else {
        pInit->alpha = atanf((points[furthest][1] - pInit->Xc[1])/(points[furthest][0] - pInit->Xc[0]));
    }
    pInit->a = max_dist;
    pInit->b = min_dist;
}

void mass_center(real **points, int N, real C[2], real *pR)
{
    int i;
    real d[2];
    
    C[0] = C[1] = 0.0f;
    for (i = 0; i < N; i++) {
        C[0] += points[i][0];
        C[1] += points[i][1];
    }
    C[0] /= N; C[1] /= N;
    
    *pR = 0;
    for (i = 0; i < N; i++) {
        vsub(C, points[i], d, 2);
        *pR += norm(d, 2);
    }
    *pR /= N;
}

int ellipse_fitting(real **points, int N, struct ellipse *pInit, struct ellipse *pRes)
{
    int i;
    int iter_count = 1;
    real Xc[2], a, b, alpha;
    real x[2], d, l;
    uint8_t iscircle = 0;
    
    real J_data[MAX_POINTS_NUM * PARAMS_COUNT], Jt_data[MAX_POINTS_NUM * PARAMS_COUNT];
    real e_data[MAX_POINTS_NUM];
    real Js_data[PARAMS_COUNT * PARAMS_COUNT], es_data[PARAMS_COUNT];
    
    matrix J   = {N,               PARAMS_COUNT,   J_data};
    matrix Jt  = {PARAMS_COUNT,    N,              Jt_data};
    matrix e   = {N,               1,              e_data};
    matrix Js  = {PARAMS_COUNT,    PARAMS_COUNT,   Js_data};
    matrix es  = {PARAMS_COUNT,    1,              es_data};
    
    uint8_t Xc_ok[2], alpha_ok;
    uint8_t rel_prec_achieved;
    
    Xc[0]   = pInit->Xc[0];
    Xc[1]   = pInit->Xc[1];
    a       = pInit->a;
    b       = pInit->b;
    alpha   = pInit->alpha;
    
    while (iter_count <= MAX_ITER_COUNT) {
        for (i = 0; i < N; i++) {
            global_to_canonical(points[i], alpha, Xc, x);
            if (i == 89) {
                i = 89;
            }
            if (!dist_to_ellipse(a, b, x, &d, &l)) {
                return 0;
            }
            e_data[i] = -d;
            
            iscircle = fabsr(1 - fabsr(a/b)) < FIT_REL_PREC;
            
            if (iscircle) {
                fill_jacobian(J_data + i*(PARAMS_COUNT-1), x, d, l, a, b, alpha);
            } else {
                fill_jacobian(J_data + i*PARAMS_COUNT, x, d, l, a, b, alpha);
            }
        }
        
        // solve linear system J * da = e
        // to find ellipse parameters' changes
        if (iscircle) {
            J.numCols   = PARAMS_COUNT - 1;
            Jt.numRows  = PARAMS_COUNT - 1;
            Js.numRows  = PARAMS_COUNT - 1; Js.numCols = PARAMS_COUNT- 1;
            es.numRows  = PARAMS_COUNT - 1;
            
            mtrans(&J, &Jt);
            mmult(&Jt, &J, &Js);
            mmult(&Jt, &e, &es);
            
            if (!linsolve(&Js, &es)) {
                return 0;
            }
        } else {
            J.numCols   = PARAMS_COUNT;
            Jt.numRows  = PARAMS_COUNT;
            Js.numRows  = PARAMS_COUNT; Js.numCols = PARAMS_COUNT;
            es.numRows  = PARAMS_COUNT;
            
            mtrans(&J, &Jt);
            mmult(&Jt, &J, &Js);
            mmult(&Jt, &e, &es);
            
            if (!linsolve(&Js, &es)) {
                return 0;
            }
            alpha += es_data[4];
        }
        // solution is in the vector es
        Xc[0]   += es_data[0];
        Xc[1]   += es_data[1];
        a       += es_data[2];
        b       += es_data[3];
        
        memset(Xc_ok, 0, 2 * sizeof(uint8_t));
        for (i = 0; i < 2; i++) {
            Xc_ok[i] = check_pot_zero_param(Xc[i], es_data[i]);
        }
        
        alpha_ok = 0;
        if (es.numRows < 5)
            alpha_ok = 1;
        else if (check_pot_zero_param(alpha, es_data[4]))
            alpha_ok = 1;
        
        rel_prec_achieved = Xc_ok[0] && Xc_ok[1] && alpha_ok &&
        fabsr(es_data[2]/a) < FIT_REL_PREC && fabsr(es_data[3]/b) < FIT_REL_PREC;
        
        if (rel_prec_achieved || norm(es_data, es.numRows) < powr(FIT_REL_PREC, 1.5)) {
            pRes->Xc[0] = Xc[0];
            pRes->Xc[1] = Xc[1];
            pRes->a = a;
            pRes->b = b;
            pRes->alpha = alpha;
            return 1;
        }
        
        iter_count++;
    }
    return 0;
}

void fill_jacobian(real *Ji, real x[2], real d, real l, real a, real b, real alpha)
{
    real dg_dl, dg_dx, dg_dy, dg_da, dg_db;
    real dl_dx, dl_dy, dl_da, dl_db;
    real dx_dXc, dx_dYc, dx_dalpha, dx_da, dx_db;
    real dy_dXc, dy_dYc, dy_dalpha, dy_da, dy_db;
    real dl_dXc, dl_dYc, dl_dalpha;
    real da2_dXc, da2_dYc, da2_da, da2_db, da2_dalpha;
    real db2_dXc, db2_dYc, db2_da, db2_db, db2_dalpha;
    
    if (l == 0 && d == 0) {
        memset(Ji, 0, PARAMS_COUNT * sizeof(real));
        return;
    }
    
    // auxiliary variables
    dg_dl = 4*l*l*l + 6*l*l * (a*a + b*a) +
    2*l * (a*a*a*a + b*b*b*b + 4*a*a*b*b - a*a*x[0]*x[0] - b*b*x[1]*x[1]) +
    2*a*a*b*b * (a*a + b*b - x[0]*x[0] - x[1]*x[1]);
    
    dg_dx = -2*x[0]*a*a * powr(b*b + l, 2);
    
    dg_dy = -2*x[1]*b*b * powr(a*a + l, 2);
    
    dg_da = l*l*l * 4*a + l*l * 2*a * (2*a*a + 4*b*b - x[0]*x[0]) +
    l * 4*a*b*b * (2*a*a + b*b - x[0]*x[0] - x[1]*x[1]) +
    2*a*b*b * (2*a*a*b*b - b*b*x[0]*x[0] - 2*a*a*x[1]*x[1]);
    
    dg_db = l*l*l * 4*b + l*l * 2*b * (2*b*b + 4*a*a - x[1]*x[1]) +
    l * 4*a*a*b * (a*a + 2*b*b - x[0]*x[0] - x[1]*x[1]) +
    2*a*a*b * (2*a*a*b*b - 2*b*b*x[0]*x[0] - a*a*x[1]*x[1]);
    
    
    dl_dx = - dg_dx / dg_dl;
    dl_dy = - dg_dy / dg_dl;
    dl_da = - dg_da / dg_dl;
    dl_db = - dg_db / dg_dl;
    
    dx_dXc = -cosr(alpha);
    dx_dYc = -sinr(alpha);
    dx_dalpha = x[1];
    dx_da = 0;
    dx_db = 0;
    
    dy_dXc = sinr(alpha);
    dy_dYc = -cosr(alpha);
    dy_dalpha = -x[0];
    dy_da = 0;
    dy_db = 0;
    
    dl_dXc = dl_dx * dx_dXc + dl_dy * dy_dXc;
    dl_dYc = dl_dx * dx_dYc + dl_dy * dy_dYc;
    dl_dalpha = dl_dx * dx_dalpha + dl_dy * dy_dalpha;
    
    da2_dXc = 0;
    da2_dYc = 0;
    da2_da = 2*a;
    da2_db = 0;
    da2_dalpha = 0;
    
    db2_dXc = 0;
    db2_dYc = 0;
    db2_da = 0;
    db2_db = 2*b;
    db2_dalpha = 0;
    
    // dd_dXc
    Ji[0] = (l/d) * (
                     x[0] / powr(a*a + l, 3) * (a*a*x[0] * dl_dXc + l*(a*a + l) * dx_dXc - l*x[0] * da2_dXc)
                     +
                     x[1] / powr(b*b + l, 3) * (b*b*x[1] * dl_dXc + l*(b*b + l) * dy_dXc - l*x[1] * db2_dXc)
                     );
    
    // dd_dYc
    Ji[1] = (l/d) * (
                     x[0] / powr(a*a + l, 3) * (a*a*x[0] * dl_dYc + l*(a*a + l) * dx_dYc - l*x[0] * da2_dYc)
                     +
                     x[1] / powr(b*b + l, 3) * (b*b*x[1] * dl_dYc + l*(b*b + l) * dy_dYc - l*x[1] * db2_dYc)
                     );
    
    // dd_da
    Ji[2] = (l/d) * (
                     x[0] / powr(a*a + l, 3) * (a*a*x[0] * dl_da + l*(a*a + l) * dx_da - l*x[0] * da2_da)
                     +
                     x[1] / powr(b*b + l, 3) * (b*b*x[1] * dl_da + l*(b*b + l) * dy_da - l*x[1] * db2_da)
                     );
    
    // dd_db
    Ji[3] = (l/d) * (
                     x[0] / powr(a*a + l, 3) * (a*a*x[0] * dl_db + l*(a*a + l) * dx_db - l*x[0] * da2_db)
                     +
                     x[1] / powr(b*b + l, 3) * (b*b*x[1] * dl_db + l*(b*b + l) * dy_db - l*x[1] * db2_db)
                     );
    
    // dd_alpha
    Ji[4] = (l/d) * (
                     x[0] / powr(a*a + l, 3) * (a*a*x[0] * dl_dalpha + l*(a*a + l) * dx_dalpha - l*x[0] * da2_dalpha)
                     +
                     x[1] / powr(b*b + l, 3) * (b*b*x[1] * dl_dalpha + l*(b*b + l) * dy_dalpha - l*x[1] * db2_dalpha)
                     );
}

void global_to_canonical(real X[2], real alpha, real Xc[2], real x[2])
{
    real ca = cosr(alpha);
    real sa = sinr(alpha);
    x[0] = ca * (X[0] - Xc[0]) + sa * (X[1] - Xc[1]);
    x[1] = -sa * (X[0] - Xc[1]) + ca * (X[1] - Xc[1]);
}

void canonical_to_global(real x[2], real alpha, real Xc[2], real X[2])
{
    real ca = cosr(alpha);
    real sa = sinr(alpha);
    X[0] = ca * x[0] + -sa * x[1] + Xc[0];
    X[1] = sa * x[0] + ca * x[1] + Xc[1];
}

// check potentially zero parameter
uint8_t check_pot_zero_param(real value, real delta)
{
    if (fabsr(value) < FIT_REL_PREC) {  // == 0
        return (fabsr(delta) < FIT_REL_PREC) ? 1 : 0;
    } else {
        return (fabsr(delta/value) < FIT_REL_PREC) ? 1 : 0;
    }
}

/*
 * LU decomposition of matrix A
 * Matrices L and U are return in the respective parts of the matrix A
 * indx - row perturbations array
 * Return 1 on success, 0 if the matrix is singular
 */
int ludcmp(real **a, int n, int *indx)
{
    int i, imax = 0, j, k;
    real big, dum, sum, temp;
    real vv[MAX_DCMP_N];
    
    // calculate the biggest elements in every row
    for (i = 0; i < n; i++) {
        big = 0.0f;
        for (j = 0; j < n; j++) {
            if ((temp = fabsr(a[i][j])) > big)
                big = temp;
        }
        if (big == 0.0f)
            return 0; // Singular matrix in routine ludcmp
        vv[i] = (1.0f / big);
    }
    
    for (j = 0; j < n; j++) {
        // calculate elements of U
        for (i = 0; i < j; i++) {
            sum = a[i][j];
            for (k = 0; k < i; k++)
                sum -= (a[i][k] * a[k][j]);
            a[i][j] = sum;
        }
        // calculate elements of L
        big = 0.0f;
        for (i = j; i < n; i++) {
            sum = a[i][j];
            for (k = 0; k < j; k++)
                sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
            
            // find the row with the biggest jth element
            if ((dum = vv[i] * fabsr(sum)) >= big) {
                big = dum;
                imax = i;
            }
        }
        // figure out if we need to interchange rows
        if (j != imax) {
            for (k = 0; k < n; k++) {
                dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dum;
            }
            vv[imax] = vv[j];
        }
        // remember where the jth row goes
        indx[j] = imax;
        
        // carry out division by the revealed biggest jj-th element
        dum = 1.0f / a[j][j];
        for (i = j + 1; i < n; i++)
            a[i][j] *= dum;
    }
    //free(vv);
    return 1;
}

/*
 * Solves the set of linear equations Ax = b
 * a     - matrix LU-decomposition returned by ludcmp
 * indx  - vector of row perturbations returned by ludcmp
 * b     - right side vector and solution is returned in it
 */
void lubksb(real **a, int n, real *b, int *indx)
{
    int i, ii = -1, ip, j;
    real sum;
    
    for (i = 0; i < n; i++) {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if (ii >= 0) {
            for (j = ii; ii < i; ii++)
                sum -= a[i][j] * b[j];
        }
        else if (sum)
            ii = i;
        b[i] = sum;
    }
    for (i = n - 1; i >= 0; i--) {
        sum = b[i];
        for (j = i + 1; j < n; j++)
            sum -= a[i][j] * b[j];
        b[i] = sum / a[i][i];
    }
}

/*
 * Solves system of linear algebraic equations Ax = b
 * Returns the solution in matrix b
 * Returns 1 on success
 */
int linsolve(matrix *pA, matrix *pb)
{
    real *a[MAX_DCMP_N], *pdata;
    int i, n, indx[MAX_DCMP_N];
    
    n = pA->numRows;
    pdata = pA->pData;
    
    for (i = 0; i < n; i++) {
        a[i] = &(pdata[i*n]);
    }
    
    if (!ludcmp(a, n, indx))
        return 0;
    
    lubksb(a, n, pb->pData, indx);
    return 1;
}