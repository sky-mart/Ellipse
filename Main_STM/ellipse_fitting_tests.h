//
//  ellipse_fitting_tests.h
//  Eq4SinglePrec
//
//  Created by ???? ??????? on 28.02.17.
//  Copyright © 2017 Mocsmart. All rights reserved.
//

#ifndef ellipse_fitting_tests_h
#define ellipse_fitting_tests_h

#include "settings.h"

int test_solve2(real a, real b, real c);
int testset_solve2(void);
int test_solve3(real a, real b, real c, real d);
int testset_solve3(void);
int test_solve4(real a, real b, real c, real d, real e);
int testset_solve4(void);
int testset_solve(void);

int test_dist_to_ellipse(real a, real b, real x[2], real trued);
int testset_dist_to_ellipse(void);

int test_ludcmp(mreal **a, int n);
int testset_ludcmp(void);
int test_linsolve(matrix *pA, matrix *pb);
int testset_linsolve(void);

int test_ellipse_fitting(int points_num, real a, real alpha, real noise_level);
int testset_ellipse_fitting(void);

#endif /* ellipse_fitting_tests_h */
