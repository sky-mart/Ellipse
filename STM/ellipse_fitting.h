#ifndef ELLIPSE_FITTING_H
#define ELLIPSE_FITTING_H

#include "complex.h"

typedef float complex complexf;

void cparts(complexf z, float *pr, float *pi);

int solve2(float a, float b, float c, complexf roots[2]);
int solve3(float a, float b, float c, float d, complexf roots[3]);
int solve4(float a, float b, float c, float d, float e, complexf roots[4]);

int dist_to_ellipse(float a, float b, float x[2], float *pd, float *pl);

#endif // ELLIPSE_FITTING
