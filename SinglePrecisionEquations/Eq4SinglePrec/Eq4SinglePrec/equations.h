//
//  equations.h
//  Eq4SinglePrec
//
//  Created by Влад Агиевич on 27.02.17.
//  Copyright © 2017 Mocsmart. All rights reserved.
//

#ifndef equations_h
#define equations_h

#include <complex.h>    // for complex numbers

//#define SINGLE_PREC

#ifdef SINGLE_PREC
    typedef float real;
    typedef float complex complexr;

    #define crealr  crealf
    #define cimagr  cimagf
    #define sqrtr   sqrtf
    #define powr    powf
    #define cbrtr   cbrtf
    #define cpowr   cpowf
    #define fminr   fminf
    #define fabsr   fabsf
#else
    typedef double real;
    typedef double complex complexr;

    #define crealr  creal
    #define cimagr  cimag
    #define sqrtr   sqrt
    #define powr    pow
    #define cbrtr   cbrt
    #define cpowr   cpow
    #define fminr   fmin
    #define fabsr   fabs
#endif

void cparts(complexr z, real *pr, real *pi);

int solve2(real a, real b, real c, complexr roots[2]);
int solve3(real a, real b, real c, real d, complexr roots[3]);
int solve4(real a, real b, real c, real d, real e, complexr roots[4]);

int dist_to_ellipse(real a, real b, real x[2], real *pd, real *pl);

#endif /* equations_h */
