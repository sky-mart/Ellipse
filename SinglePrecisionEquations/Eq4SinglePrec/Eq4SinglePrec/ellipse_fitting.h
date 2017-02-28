//
//  ellipse_fitting.h
//  Eq4SinglePrec
//
//  Created by Влад Агиевич on 28.02.17.
//  Copyright © 2017 Mocsmart. All rights reserved.
//

#ifndef ellipse_fitting_h
#define ellipse_fitting_h

#include <complex.h>    // for complex numbers

//#define SINGLE_PREC_DEBUG
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
    #define sinr    sinf
    #define cosr    cosf
    #define cabsr   cabsf
    #define hypotr  hypotf
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
    #define sinr    sin
    #define cosr    cos
    #define cabsr   cabs
    #define hypotr  hypot
#endif

struct ellipse
{
    real Xc[2];
    real a;
    real b;
    real alpha;
};

// for possible arm_math implementation
typedef struct
{
    uint16_t numRows;
    uint16_t numCols;
    real *pData;
} mat;

typedef mat matrix;
#define dotprod     dotprodr
#define vsub        vec_sub
#define vmax        vec_max
#define mtrans      mat_trans
#define mmult       mat_mult
#define minit       mat_init
#define msub        mat_sub
////////////////////////////////////////

// algorithm and tests parameters
#define MAX_POINTS_NUM              500
#define MAX_DCMP_N                  5
#define MAX_ITER_COUNT              50
#define MAX_ERROR_INCREASE_COUNT    4
#define PARAMS_COUNT                5
#define FIT_REL_PREC                1e-4f

void cparts(complexr z, real *pr, real *pi);

int solve2(real a, real b, real c, complexr roots[2]);
int solve3(real a, real b, real c, real d, complexr roots[3]);
int solve4(real a, real b, real c, real d, real e, complexr roots[4]);

int dist_to_ellipse(real a, real b, real x[2], real *pd, real *pl);

void mass_center(real **points, int N, real C[2], real *pR);
void fill_jacobian(real *Ji, real x[2], real d, real l, real a, real b, real alpha);
void global_to_canonical(real X[2], real alpha, real Xc[2], real x[2]);
void canonical_to_global(real x[2], real alpha, real Xc[2], real X[2]);
void ellipse_fitting_init_guess(real **points, int N, struct ellipse *pInit);
int ellipse_fitting(real **points, int N, struct ellipse *pInit, struct ellipse *pRes);
uint8_t check_pot_zero_param(real value, real delta);
int linsolve(matrix *pA, matrix *pb);
int ludcmp(real **a, int n, int *indx);
void lubksb(real **a, int n, real *b, int *indx);

void vec_sub(real *a, real *b, real *d, uint n);
void vec_max(real *v, uint n, real *pResult, uint *pIndex);
void mat_init(matrix *pA, uint16_t nRows, uint16_t nCols, real *pData);
int mat_trans(matrix *pA, matrix *pAt);
int mat_sub(matrix *pA, matrix *pB, matrix *pC);
int mat_mult(matrix *pA, matrix *pB, matrix *pC);

#endif /* ellipse_fitting_h */
