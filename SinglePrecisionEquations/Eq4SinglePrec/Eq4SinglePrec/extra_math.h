//
//  extra_math.h
//  Eq4SinglePrec
//
//  Created by Влад Агиевич on 28.02.17.
//  Copyright © 2017 Mocsmart. All rights reserved.
//

#ifndef extra_math_h
#define extra_math_h

#include "settings.h"

void cparts(complexr z, real *pr, real *pi);
void cprint(complexr z);
//real mypowr(real b, real e);

void dotprodr(real *a, real *b, uint n, real *pResult);
void vec_sub(real *a, real *b, real *d, uint n);
void vec_max(real *v, uint n, real *pResult, uint *pIndex);
real norm(real *v, uint n);
real dist(real *u, real *v, int n);

void mat_init(matrix *pA, uint16_t nRows, uint16_t nCols, real *pData);
int mat_trans(matrix *pA, matrix *pAt);
int mat_sub(matrix *pA, matrix *pB, matrix *pC);
int mat_mult(matrix *pA, matrix *pB, matrix *pC);

int ludcmp(real **a, int n, int *indx);
void lubksb(real **a, int n, real *b, int *indx);
int linsolve(matrix *pA, matrix *pb);



#endif /* extra_math_h */
