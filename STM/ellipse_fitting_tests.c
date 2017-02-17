#include "ellipse_fitting.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"

#define EQ_PREC         5e-5f   // precision for equations    
#define DST_PREC        1e-6f   // precision for distance
#define MAT_PREC        1e-6f   // precision for matrices
#define ANGLE_PREC      1e-3f


float **alloc_array2d(size_t nRows, size_t nCols)
{
    float *t1, **t2;
    int i;
    if ( !(t1 = (float *)malloc(nRows * nCols * sizeof(float))) )
        return NULL;
    if ( !(t2 = (float **)malloc(nRows * sizeof(float*))) ) {
        free(t1);
        return NULL;
    }
    for (i = 0; i < nRows; i++) {
        t2[i] = &t1[i * nCols];
    }
    return t2;
}

void free_array2d(float **a)
{
    free(a[0]);
    free(a);
}

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

int test_ludcmp(float **a, int n)
{
    float A_data[MAX_DCMP_N*MAX_DCMP_N];
    float L_data[MAX_DCMP_N*MAX_DCMP_N];
    float U_data[MAX_DCMP_N*MAX_DCMP_N];
    float LxU_data[MAX_DCMP_N*MAX_DCMP_N], dum;
    arm_matrix_instance_f32 A, L, U, LxU;
    int i, j, indx[MAX_DCMP_N];
    uint32_t ii;   
    
    for (i = 0; i < n; i++) 
        for (j = 0; j < n; j++) 
            A_data[i*n + j] = a[i][j];
    
    arm_mat_init_f32(&A, n, n, A_data);
    arm_mat_init_f32(&L, n, n, L_data);
    arm_mat_init_f32(&U, n, n, U_data);
    arm_mat_init_f32(&LxU, n, n, LxU_data);
    
    if (!ludcmp(a, n, indx)) {
        return 0;
    }
    
    memset(L_data, 0, n * n * sizeof(float)); 
    memset(U_data, 0, n * n * sizeof(float)); 
    for (i = 0; i < n; i++) {
        L_data[i*n + i] = 1.0f;
        for (j = 0; j < i; j++)  
            L_data[i*n + j] = a[i][j];
        for (j = i; j < n; j++)
            U_data[i*n + j] = a[i][j];
    }

    for (i = 0; i < n; i++) {
        for (j  = 0; j < n; j++) {
            dum = A_data[i*n + j];
            A_data[i*n + j] = A_data[indx[i]*n + j];
            A_data[indx[i]*n + j] = dum;
        }
    }
    
    arm_mat_mult_f32(&L, &U, &LxU);
    arm_mat_sub_f32(&LxU, &A, &L);
    arm_max_f32(L_data, n * n, &dum, &ii);
    
    return fabs(dum) < MAT_PREC;
}

int testset_ludcmp()
{
    int result = 1;
    float **a;
    
    if ( !(a = alloc_array2d(5, 5)) ) {
        return 0;
    }
    a[0][0] = 1.0f;     a[0][1] = 2.0f;     a[0][2] = 5.0f;
    a[1][0] = 0.0f;     a[1][1] = 3.0f;     a[1][2] = 21.0f;
    a[2][0] = 7.0f;     a[2][1] = 0.0f;     a[2][2] = 2.0f;
    result &= test_ludcmp(a, 3);
    
    a[0][0] = -1.0f;    a[0][1] = 2.0f;     a[0][2] = 15.0f;    a[0][3] = 2.5f;
    a[1][0] = 3.0f;     a[1][1] = 4.0f;     a[1][2] = -21.0f;   a[1][3] = -1.3f;
    a[2][0] = -7.0f;    a[2][1] = 0.0f;     a[2][2] = 1.0f;     a[2][3] = 0.0f;
    a[3][0] = -2.0f;    a[3][1] = 3.0f;     a[3][2] = 5.0f;     a[3][3] = 4.0f;
    result &= test_ludcmp(a, 4);
    
    a[0][0] = 1.5f;     a[0][1] = 2.0f;     a[0][2] = -0.4f;    a[0][3] = 2.5f;     a[0][4] = 8.8f;
    a[1][0] = 3.0f;     a[1][1] = 4.0f;     a[1][2] = -21.0f;   a[1][3] = 100.0f;   a[1][4] = 0.1f;
    a[2][0] = -7.1f;    a[2][1] = 0.9f;     a[2][2] = 1.0f;     a[2][3] = 0.0f;     a[2][4] = 1.0f;
    a[3][0] = -2.2f;    a[3][1] = 3.0f;     a[3][2] = 5.5f;     a[3][3] = 4.7f;     a[3][4] = 72.0f;
    a[4][0] = 0.3f;     a[4][1] = 3.8f;     a[4][2] = 7.4f;     a[4][3] = 4.7f;     a[4][4] = -26.0f;
    result &= test_ludcmp(a, 5);
    
    free_array2d(a);
    return result;
}

int test_linsolve(arm_matrix_instance_f32 *pA, arm_matrix_instance_f32 *pb)
{
    float A_src_data[MAX_DCMP_N*MAX_DCMP_N];
    float b_src_data[MAX_DCMP_N], b_acq_data[MAX_DCMP_N];
    arm_matrix_instance_f32 A_src, b_src, b_acq;
    int i, n;
    float dum, big;
    
    n = pA->numRows;
    
    for (i = 0; i < n; i++) {
        A_src_data[i] = pA->pData[i];
        b_src_data[i] = pb->pData[i];
    }
    for (i = n; i < n * n; i++) {
        A_src_data[i] = pA->pData[i];
    }
    arm_mat_init_f32(&A_src, n, n, A_src_data);
    arm_mat_init_f32(&b_src, n, 1, b_src_data);
    arm_mat_init_f32(&b_acq, n, 1, b_acq_data);
    
    if (!linsolve(pA, pb))
        return 0;

    arm_mat_mult_f32(&A_src, pb, &b_acq);
    arm_mat_sub_f32(&b_src, &b_acq, pb);
    
    big = 0.0f;
    for (i = 0; i < n; i++) {
        if ( (dum = fabsf(pb->pData[i])) > big) {
            big = dum;
        }
    }
    return fabsf(big) < MAT_PREC;
}

int testset_linsolve()
{
    float A_data[MAX_DCMP_N*MAX_DCMP_N];
    float b_data[MAX_DCMP_N];
    arm_matrix_instance_f32 A, b;
    int n = 3;
    int result = 1;
    
    A_data[0*n + 0] = 1.0f;     A_data[0*n + 1] = 2.0f;     A_data[0*n + 2] = 5.0f;
    A_data[1*n + 0] = 0.0f;     A_data[1*n + 1] = 3.0f;     A_data[1*n + 2] = 21.0f;
    A_data[2*n + 0] = 7.0f;     A_data[2*n + 1] = 0.0f;     A_data[2*n + 2] = 2.0f;
    
    b_data[0] = 3.0f;           b_data[1] = 1.0f;           b_data[2] = 2.5f;
    
    arm_mat_init_f32(&A, n, n, A_data);
    arm_mat_init_f32(&b, n, 1, b_data);
    
    result &= testset_ludcmp();
    result &= test_linsolve(&A, &b);
    return result;
}

float randf(float from, float to)
{
    float val = ((float)rand()) / RAND_MAX;
    return (to - from) * val + from;
}

void generate_points(int points_num, struct ellipse *pEl, float noise_level, float **points)
{
    float step = 2 * M_PI / points_num;
    float t, x[2];
    int i;

    t = 0;
    for (i = 0; i < points_num; i++) {
        x[0] = pEl->a * cosf(t) * (1 + randf(-noise_level, noise_level));
        x[1] = pEl->b * sinf(t) * (1 + randf(-noise_level, noise_level));
        canonical_to_global(x, pEl->alpha, pEl->Xc, points[i]);
        t += step;
    }
}

float angle_to_0_2pi(float x)
{
    while (x < 0) 
        x += 2 * M_PI;
    while (x >= 2 * M_PI)
        x -= 2 * M_PI;
    return x;
}

int angles_cmp(float alpha, float beta)
{
    if (fabsf(alpha) < ANGLE_PREC)
        alpha = 0;
    if (fabsf(beta) < ANGLE_PREC)
        beta = 0;
    alpha   = angle_to_0_2pi(alpha);
    beta    = angle_to_0_2pi(beta);
    return fabsf(alpha - beta) < ANGLE_PREC;
}

int ellipse_cmp(struct ellipse *pSrc, struct ellipse *pRes)
{
    int i;
    uint8_t aeqa, beqb; // src_param is equal res_param
    uint8_t aeqb, beqa;
    uint8_t srciscircle;
    for (i = 0; i < 2; i++) {
        if (pSrc->Xc[i] == 0) {
            if (fabsf(pRes->Xc[i]) > FIT_REL_PREC)
                return 0;
        } else if (fabsf(1 - pRes->Xc[i] / pSrc->Xc[i]) > FIT_REL_PREC)
            return 0;
    }
    
    aeqa = fabsf(1.0f - pRes->a / pSrc->a) < FIT_REL_PREC;
    beqb = fabsf(1.0f - pRes->b / pSrc->b) < FIT_REL_PREC;
    srciscircle = fabsf(1.0f - pSrc->a / pSrc->b) < FIT_REL_PREC;
    if (aeqa && beqb) {
        if (!srciscircle)
            return angles_cmp(pRes->alpha, pSrc->alpha) || 
                angles_cmp(pRes->alpha, pSrc->alpha + M_PI);
        return 1;
    } else {
        aeqb = fabsf(1.0f - pRes->b / pSrc->a) < FIT_REL_PREC;
        beqb = fabsf(1.0f - pRes->a / pSrc->b) < FIT_REL_PREC;
        if (aeqb && beqa) {
            if (srciscircle)
                return angles_cmp(pRes->alpha, pSrc->alpha + M_PI / 2.0f) || 
                    angles_cmp(pRes->alpha, pSrc->alpha + 3.0f * M_PI / 2.0f);
            return 1;
        }
    }
    return 0;
}

int test_ellipse_fitting(int points_num, float a, float alpha, float noise_level)
{
    float **points;
    struct ellipse init, src, res;
 
    src.Xc[0] = 1.0;
    src.Xc[1] = 1.0;
    src.a = a;
    src.b = 1.0;
    src.alpha = alpha;
    
    if ( !(points = alloc_array2d(points_num, 2)) )
        return 0;
    
    generate_points(points_num, &src, noise_level, points);
    ellipse_fitting_init_guess(points, points_num, &init);
    if (!ellipse_fitting(points, points_num, &init, &res)) {
        free_array2d(points);
        return 0;
    }
    
    free_array2d(points);
    return ellipse_cmp(&src, &res);
}

int testset_ellipse_fitting()
{
    int result = 1;
    result &= test_ellipse_fitting(500, 1.5f, M_PI / 6.0f, 0);
    return result;
}
