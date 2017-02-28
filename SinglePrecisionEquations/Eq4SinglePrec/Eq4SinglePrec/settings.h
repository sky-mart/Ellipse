//
//  settings.h
//  Eq4SinglePrec
//
//  Created by Влад Агиевич on 28.02.17.
//  Copyright © 2017 Mocsmart. All rights reserved.
//

#ifndef settings_h
#define settings_h

#include <complex.h>    
#include <math.h>

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

// algorithm and tests parameters
#define MAX_POINTS_NUM              500
#define MAX_DCMP_N                  5
#define MAX_ITER_COUNT              50
#define MAX_ERROR_INCREASE_COUNT    4
#define PARAMS_COUNT                5
#define FIT_REL_PREC                1e-4f


#endif /* settings_h */
