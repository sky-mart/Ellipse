//
//  tests.c
//  Equation4_x86
//
//  Created by Влад Агиевич on 06.03.17.
//  Copyright © 2017 Mocsmart. All rights reserved.
//

#include "tests.h"

void cparts(rcomplex z, real *pr, real *pi)
{
    *pr = crealr(z);
    *pi = crealr(z);
}

int test_solve2(real a, real b, real c,
                int (*solve2)(real, real, real, rcomplex[2]))
{
    rcomplex roots[2], z;
    real rp[3], ip[3];
    int i;
    if ((*solve2)(a, b, c, roots)) {
        for (i = 0; i < 2; i++) {
            z = roots[i];
            cparts(z, &rp[i], &ip[i]);
            if (isnan(rp[i]) || isnan(ip[i]) ||
                cabsr(a*z*z + b*z + c) >= EQ_PREC) {
                return 0;
            }
        }
    }
    return 1;
}

int testset_solve2(int (*solve2)(real, real, real, rcomplex[2]))
{
    int result = 1;
    result &= test_solve2(1.0f, 5.0f, 6.0f, solve2);
    result &= test_solve2(1.0f, 1.0f, 2.5f, solve2);
    result &= test_solve2(1.0f, -4.0f, 4.0f, solve2);
    result &= test_solve2(0.0f, 5.0f, 5.0f, solve2);
    result &= test_solve2(0.0f, 0.0f, 3.0f, solve2);
    result &= test_solve2(0.0f, 3.0f, 0.0f, solve2);
    result &= test_solve2(1.0f, 5.0f, 0.0f, solve2);
    result &= test_solve2(1.0f, 0.0f, 9.0f, solve2);
    return result;
}

int test_solve3(real a, real b, real c, real d,
                int (*solve3)(real, real, real, real, rcomplex[3]))
{
    rcomplex roots[3], z;
    real rp[3], ip[3];
    int i;
    if ((*solve3)(a, b, c, d, roots)) {
        for (i = 0; i < 3; i++) {
            z = roots[i];
            cparts(z, &rp[i], &ip[i]);
            if (isnan(rp[i]) || isnan(ip[i]) ||
                cabsr(a*z*z*z + b*z*z + c*z + d) >= EQ_PREC) {
                return 0;
            }
        }
    }
    return 1;
}

int testset_solve3(int (*solve3)(real, real, real, real, rcomplex[3]))
{
    int result = 1;
    result &= test_solve3(1.0f, 5.0f, 6.0f, -3.0f, solve3);
    result &= test_solve3(1.0f, 3.0f, -3.0f, -1.0f, solve3);
    result &= test_solve3(1.0f, 6.0f, 9.0f, -1.0f, solve3);
    result &= test_solve3(1.0f, 6.0f, 9.0f, 22.0f, solve3);
    return result;
}

int test_solve4(real a, real b, real c, real d, real e,
                int (*solve4)(real, real, real, real, real, rcomplex[4]))
{
    rcomplex roots[4], z;
    real rp[4], ip[4];
    int i;
    if ((*solve4)(a, b, c, d, e, roots)) {
        for (i = 0; i < 4; i++) {
            z = roots[i];
            cparts(z, &rp[i], &ip[i]);
            if (isnan(rp[i]) || isnan(ip[i]) ||
                cabsr(a*z*z*z*z + b*z*z*z + c*z*z + d*z + e) >= EQ_PREC) {
                return 0;
            }
        }
    }
    return 1;
}

int testset_solve4(int (*solve4)(real, real, real, real, real, rcomplex[4]))
{
    int result = 1;
    result &= test_solve4(1.0f, 5.0f, 6.0f, -3.0f, -21.0f, solve4);
    return result;
}

int testset_solve(int (*solve2)(real, real, real, rcomplex[2]),
                  int (*solve3)(real, real, real, real, rcomplex[3]),
                  int (*solve4)(real, real, real, real, real, rcomplex[4]))
{
    int result = 1;
    result &= testset_solve2(solve2);
    result &= testset_solve3(solve3);
    result &= testset_solve4(solve4);
    return result;
}
