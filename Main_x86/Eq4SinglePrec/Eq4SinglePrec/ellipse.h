//
//  ellipse.h
//  Eq4SinglePrec
//
//  Created by Влад Агиевич on 28.02.17.
//  Copyright © 2017 Mocsmart. All rights reserved.
//

#ifndef ellipse_h
#define ellipse_h

#include "settings.h"

struct ellipse
{
    real Xc[2];
    real a;
    real b;
    real alpha;
};

int dist_to_ellipse(real a, real b, real x[2], real *pd, real *pl);

void ellipse_fitting_init_guess(real **points, int N, struct ellipse *pInit);
void mass_center(real **points, int N, real C[2], real *pR);

int ellipse_fitting(real **points, int N, struct ellipse *pInit, struct ellipse *pRes);
void fill_jacobian(mreal *Ji, real x[2], real d, real l, real a, real b, real alpha);
void global_to_canonical(real X[2], real alpha, real Xc[2], real x[2]);
void canonical_to_global(real x[2], real alpha, real Xc[2], real X[2]);
uint8_t check_pot_zero_param(real value, real delta);

#endif /* ellipse_h */
