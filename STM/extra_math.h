//
//  extra_math.h
//  Eq4SinglePrec
//
//  Created by ???? ??????? on 28.02.17.
//  Copyright © 2017 Mocsmart. All rights reserved.
//

#ifndef extra_math_h
#define extra_math_h

#include "settings.h"

#define SQRT3   1.732050807568877f
#define M_PI    3.1415926535897932384626433832795f

void cparts(complexr z, real *pr, real *pi);
void cprint(complexr z);
//real mypowr(real b, real e);

void dotprodf(float *a, float *b, uint n, float *pResult);
void dotprod(double *a, double *b, uint n, double *pResult);
void vec_sub(real *a, real *b, real *d, uint n);
void vec_maxf(float *v, uint n, float *pResult, uint *pIndex);
void vec_max(double *v, uint n, double *pResult, uint *pIndex);
float normf(float *v, uint n);
double norm(double *v, uint n);
real dist(real *u, real *v, int n);

void mat_init(matrix *pA, uint16_t nRows, uint16_t nCols, mreal *pData);
int mat_trans(matrix *pA, matrix *pAt);
int mat_sub(matrix *pA, matrix *pB, matrix *pC);
int mat_mult(matrix *pA, matrix *pB, matrix *pC);

int ludcmp(mreal **a, int n, int *indx);
void lubksb(mreal **a, int n, mreal *b, int *indx);
int linsolve(matrix *pA, matrix *pb);



#endif /* extra_math_h */
