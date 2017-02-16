#ifndef ELLIPSE_FITTING_H
#define ELLIPSE_FITTING_H

#include "complex.h"    // for complex numbers
#include "arm_math.h"   // for matrix operations

typedef float complex complexf;

void cparts(complexf z, float *pr, float *pi);

int solve2(float a, float b, float c, complexf roots[2]);
int solve3(float a, float b, float c, float d, complexf roots[3]);
int solve4(float a, float b, float c, float d, float e, complexf roots[4]);

int dist_to_ellipse(float a, float b, float x[2], float *pd, float *pl);

#define MAX_DCMP_N 5
int ludcmp(float **a, int n, int *indx);
void lubksb(float **a, int n, float *b, int *indx);
int linsolve(arm_matrix_instance_f32 *pA, arm_matrix_instance_f32 *pb);

#endif // ELLIPSE_FITTING
