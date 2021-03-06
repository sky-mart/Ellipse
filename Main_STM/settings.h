//
//  settings.h
//  Eq4SinglePrec
//
//  Created by ???? ??????? on 28.02.17.
//  Copyright � 2017 Mocsmart. All rights reserved.
//

#ifndef settings_h
#define settings_h

#include <complex.h>    
#include <math.h>
#include <arm_math.h>
#include <stdint.h>

typedef unsigned int uint;

//#define SINGLE_PREC_DEBUG
#define SINGLE_PREC
#define MAT_SINGLE_PREC

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
    #define normr   normf
    #define alloc_array2dr  alloc_array2df
    #define free_array2dr   free_array2df
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
    #define normr   norm
    #define alloc_array2dr  alloc_array2d
    #define free_array2dr   free_array2d
#endif

#ifdef MAT_SINGLE_PREC
    typedef float mreal;
    #define dotprodm        mdotprodf
    #define normm           normf
    #define vmaxm           vec_maxf
    #define alloc_array2dm  alloc_array2df
    #define free_array2dm   free_array2df
#define fabsm   fabsf
#else
    typedef double mreal;
    #define dotprodm        mdotprod
    #define normm           norm
    #define vmaxm           vec_max
    #define alloc_array2dm  alloc_array2d
    #define free_array2dm   free_array2d
#define fabsm   fabs
#endif


#ifndef ARM_MATRICES
    typedef struct
    {
        uint16_t numRows;
        uint16_t numCols;
        mreal *pData;
    } mat;

    typedef mat matrix;

    #define vsub        vec_sub
    #define mtrans      mat_trans
    #define mmult       mat_mult
    #define minit       mat_init
    #define msub        mat_sub
#else
    typedef arm_matrix_instance_f32 matrix;
    #define vsub        arm_sub_f32
    #define mtrans      arm_mat_trans_f32
    #define mmult       arm_mat_mult_f32
    #define minit       arm_mat_init_f32
    #define msub        arm_mat_sub_f32
#endif

    
// algorithm and tests parameters
#define MAX_POINTS_NUM              500
#define MAX_DCMP_N                  5
#define MAX_ITER_COUNT              50
#define MAX_ERROR_INCREASE_COUNT    4
#define PARAMS_COUNT                5
#define FIT_REL_PREC                1e-3f
#define FIT_ABS_PREC                1e-2f


#endif /* settings_h */
