#ifndef ELLIPSE_FITTING_TESTS_H
#define ELLIPSE_FITTING_TESTS_H

#include "arm_math.h"

int test_solve2(float a, float b, float c);
int testset_solve2(void);
int test_solve3(float a, float b, float c, float d);
int testset_solve3(void);
int test_solve4(float a, float b, float c, float d, float e);
int testset_solve4(void);
int testset_solve(void);

int test_dist_to_ellipse(float a, float b, float x[2], float trued);
int testset_dist_to_ellipse(void);

int test_ludcmp(float **a, int n);
int testset_ludcmp(void);
int test_linsolve(arm_matrix_instance_f32 *pA, arm_matrix_instance_f32 *pb);
int testset_linsolve(void);

#endif // ELLIPSE_FITTING_TESTS_H
