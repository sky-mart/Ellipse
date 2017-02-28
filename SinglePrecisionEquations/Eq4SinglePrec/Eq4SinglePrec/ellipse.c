//
//  ellipse.c
//  Eq4SinglePrec
//
//  Created by Влад Агиевич on 28.02.17.
//  Copyright © 2017 Mocsmart. All rights reserved.
//

#include "ellipse.h"
#include "equations.h"
#include "extra_math.h"
#include <float.h>
#include <string.h>

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

/*
 * Initial guess for ellipse_fitting routine
 * Returns Xc, a, b, alpha in pointers
 */
void ellipse_fitting_init_guess(real **points, int N, struct ellipse *pInit)
{
    real R;
    int closest, furthest, i;
    real min_dist = 0.0f, max_dist = 0.0f, d;
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

