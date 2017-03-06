//
//  settings.h
//  Equation4_x86
//
//  Created by Влад Агиевич on 06.03.17.
//  Copyright © 2017 Mocsmart. All rights reserved.
//

#ifndef settings_h
#define settings_h

#include <complex.h>
#include <math.h>

#define SINGLE_PREC
#define DIRECT_SINGLE

#ifdef SINGLE_PREC
    typedef float real;
    typedef float complex rcomplex;
    #define crealr  crealf
    #define cimagr  cimagf
    #define sqrtr   sqrtf
    #define powr    powf
    #define cbrtr   cbrtf
    #define cpowr   cpowf
    #define fminr   fminf
    #define fabsr   fabsf
    #define acosr   acosf
    #define sinr    sinf
    #define cosr    cosf
    #define cabsr   cabsf
    #define hypotr  hypotf
    #define normr   normf
    #define alloc_array2dr  alloc_array2df
    #define free_array2dr   free_array2df
#else
    typedef double real;
    typedef double complex rcomplex;
    #define crealr  creal
    #define cimagr  cimag
    #define sqrtr   sqrt
    #define powr    pow
    #define cbrtr   cbrt
    #define cpowr   cpow
    #define fminr   fmin
    #define fabsr   fabs
    #define acosr   acos
    #define sinr    sin
    #define cosr    cos
    #define cabsr   cabs
    #define hypotr  hypot
    #define normr   norm
    #define alloc_array2dr  alloc_array2d
    #define free_array2dr   free_array2d
#endif

#ifdef DIRECT_SINGLE
    typedef float dreal;
    typedef float complex dcomplex;
    #define creald  crealf
    #define cimagd  cimagf
    #define sqrtd   sqrtf
    #define powd    powf
    #define cbrtd   cbrtf
    #define cpowd   cpowf
    #define fmind   fminf
    #define fabsd   fabsf
    #define acosd   acosf
    #define sind    sinf
    #define cosd    cosf
    #define cabsd   cabsf
    #define hypotd  hypotf
    #define normd   normf
#else
    typedef double dreal;
    typedef double complex dcomplex;
    #define creald  creal
    #define cimagd  cimag
    #define sqrtd   sqrt
    #define powd    pow
    #define cbrtd   cbrt
    #define cpowd   cpow
    #define fmind   fmin
    #define fabsd   fabs
    #define acosd   acos
    #define sind    sin
    #define cosd    cos
    #define cabsd   cabs
    #define hypotd  hypot
    #define normd   norm
#endif

#endif /* settings_h */
