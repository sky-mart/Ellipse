#include "ellipse_fitting.h"
#include "math.h"

#define EQ_PREC     5e-5f
#define DST_PREC    1e-6f

int test_solve2(float a, float b, float c)
{
    complexf roots[2], z;
    float rp[3], ip[3];
    int i;
    if (solve2(a, b, c, roots)) {
        for (i = 0; i < 2; i++) {
            z = roots[i];
            cparts(z, &rp[i], &ip[i]);
            if (isnan(rp[i]) || isnan(ip[i]) ||
                cabsf(a*z*z + b*z + c) >= EQ_PREC) {
                return 0;
            }
        }            
    }
    return 1;
}

int testset_solve2()
{
    int result = 1;
    result &= test_solve2(1.0f, 5.0f, 6.0f);
    result &= test_solve2(1.0f, 1.0f, 2.5f);
    result &= test_solve2(1.0f, -4.0f, 4.0f);
    result &= test_solve2(0.0f, 5.0f, 5.0f);
    result &= test_solve2(0.0f, 0.0f, 3.0f);
    result &= test_solve2(0.0f, 3.0f, 0.0f);
    result &= test_solve2(1.0f, 5.0f, 0.0f);
    result &= test_solve2(1.0f, 0.0f, 9.0f);
    return result;
}

int test_solve3(float a, float b, float c, float d)
{
    complexf roots[3], z;
    float rp[3], ip[3];
    int i;
    if (solve3(a, b, c, d, roots)) {
        for (i = 0; i < 3; i++) {
            z = roots[i];
            cparts(z, &rp[i], &ip[i]);
            if (isnan(rp[i]) || isnan(ip[i]) ||
                cabsf(a*z*z*z + b*z*z + c*z + d) >= EQ_PREC) {
                return 0;
            }
        }            
    }
    return 1;
}

int testset_solve3()
{
    int result = 1;
    result &= test_solve3(1.0f, 5.0f, 6.0f, -3.0f);
    result &= test_solve3(1.0f, 3.0f, -3.0f, -1.0f);
    result &= test_solve3(1.0f, 6.0f, 9.0f, -1.0f);
    result &= test_solve3(1.0f, 6.0f, 9.0f, 22.0f);
    return result;
}

int test_solve4(float a, float b, float c, float d, float e)
{
    complexf roots[4], z;
    float rp[4], ip[4];
    int i;
    if (solve4(a, b, c, d, e, roots)) {
        for (i = 0; i < 4; i++) {
            z = roots[i];
            cparts(z, &rp[i], &ip[i]);
            if (isnan(rp[i]) || isnan(ip[i]) ||
                cabsf(a*z*z*z*z + b*z*z*z + c*z*z + d*z + e) >= EQ_PREC) {
                return 0;
            }
        }            
    }
    return 1;
}

int testset_solve4()
{
    int result = 1;
    result &= test_solve4(1.0f, 5.0f, 6.0f, -3.0f, -21.0f);
    return result;
}

int testset_solve()
{
    int result = 1;
    result &= testset_solve2();
    result &= testset_solve3();
    result &= testset_solve4();
    return result;
}

int test_dist_to_ellipse(float a, float b, float x[2], float trued)
{
    float d, l;
    if (!dist_to_ellipse(a, b, x, &d, &l)) {
        return 0;
    }
    return fabsf(d - trued) < DST_PREC;
}

int testset_dist_to_ellipse()
{
    int result = 1;
    float x1[2] = {4.0f, 0.0f};
    float x2[2] = {0.0f, 5.0f};
    float x3[2] = {-3.0f, 0.0f};
    float x4[2] = {0.0f, -8.0f};
    float x5[2] = {2.0f, 2.0f};
    float x6[2] = {-3.0f, 4.0f};
    result &= test_dist_to_ellipse(2.0f, 1.0f, x1, 2.0f);
    result &= test_dist_to_ellipse(2.0f, 1.0f, x2, 4.0f);
    result &= test_dist_to_ellipse(2.0f, 1.0f, x3, 1.0f);
    result &= test_dist_to_ellipse(2.0f, 1.0f, x4, 7.0f);
    result &= test_dist_to_ellipse(1.0f, 1.0f, x5, hypotf(x5[0], x5[1]) - 1);
    result &= test_dist_to_ellipse(1.0f, 1.0f, x6, hypotf(x6[0], x6[1]) - 1);
    return result;
}
