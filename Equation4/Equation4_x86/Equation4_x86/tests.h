//
//  tests.h
//  Equation4_x86
//
//  Created by Влад Агиевич on 06.03.17.
//  Copyright © 2017 Mocsmart. All rights reserved.
//

#ifndef tests_h
#define tests_h

#include "settings.h"

#define EQ_PREC         5e-5f   // precision for equations

int test_solve2(real a, real b, real c,
                int (*solve2)(real, real, real, rcomplex[2]));

int testset_solve2(int (*solve2)(real, real, real, rcomplex[2]));

int test_solve3(real a, real b, real c, real d,
                int (*solve3)(real, real, real, real, rcomplex[3]));

int testset_solve3(int (*solve3)(real, real, real, real, rcomplex[3]));

int test_solve4(real a, real b, real c, real d, real e,
                int (*solve4)(real, real, real, real, real, rcomplex[4]));

int testset_solve4(int (*solve4)(real, real, real, real, real, rcomplex[4]));

int testset_solve(int (*solve2)(real, real, real, rcomplex[2]),
                  int (*solve3)(real, real, real, real, rcomplex[3]),
                  int (*solve4)(real, real, real, real, real, rcomplex[4]));

#endif /* tests_h */
