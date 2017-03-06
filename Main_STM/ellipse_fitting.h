#ifndef ELLIPSE_FITTING_H
#define ELLIPSE_FITTING_H

#include "complex.h"    // for complex numbers
#include "arm_math.h"   // for matrix operations

#define M_PI       3.14159265358979323846f  /* pi */

#define MAX_POINTS_NUM  500
#define MAX_DCMP_N      5
#define MAX_ITER_COUNT              50
#define MAX_ERROR_INCREASE_COUNT    4
#define PARAMS_COUNT                5
#define FIT_REL_PREC                1e-4f

typedef float complex complexf;

struct ellipse
{
    float Xc[2];
    float a;
    float b;
    float alpha;
};

void cparts(complexf z, float *pr, float *pi);

int solve2(float a, float b, float c, complexf roots[2]);
int solve3(float a, float b, float c, float d, complexf roots[3]);
int solve4(float a, float b, float c, float d, float e, complexf roots[4]);

int dist_to_ellipse(float a, float b, float x[2], float *pd, float *pl);

void canonical_to_global(float x[2], float alpha, float Xc[2], float X[2]);

int ludcmp(float **a, int n, int *indx);
void lubksb(float **a, int n, float *b, int *indx);
int linsolve(arm_matrix_instance_f32 *pA, arm_matrix_instance_f32 *pb);

void ellipse_fitting_init_guess(float **points, int N, struct ellipse *pInit);
int ellipse_fitting(float **points, int N, struct ellipse *pInit, struct ellipse *pRes);

#endif // ELLIPSE_FITTING
