#include "ellipse_fitting.h"

#include "math.h"
#include "float.h"      // FLT_MAX
#include "stdlib.h"     // malloc
#include "string.h"     // memset

void cparts(complexf z, float *pr, float *pi)
{
    *pr = crealf(z);
    *pi = cimagf(z);
}

//float mypowf(float b, float e) 
//{
//    return (b >= 0) ? powf(b, e) : -powf(-b, e);
//}

/*
* solves equation ax^2 + bx + c = 0
* puts roots in the array
* returns 1 on success
*/
int solve2(float a, float b, float c, complexf roots[2])
{
    float D;
    if (a == 0) {
        if (b == 0) {
            return 0; // any x is a solution
        }
        roots[0] = roots[1] = -c / b;
        return 1;
    }
    
    D = b*b - 4*a*c;
    if (D >= 0) {
        roots[0] = (-b + sqrtf(D)) / (2*a);
        roots[1] = (-b - sqrtf(D)) / (2*a);
    } else {
        roots[0] = (-b + sqrtf(-D) * I) / (2*a);
        roots[1] = (-b - sqrtf(-D) * I) / (2*a);
    }
    return 1;
}

/*
* solves equation ax^3 + bx^2 + cx + d = 0
* puts roots in the array
* returns 1 on success
*/
int solve3(float a, float b, float c, float d, complexf roots[3])
{
    float p, q, Q;
    complexf alpha, beta;
    // float rp[5], ip[5]; // debug
    
    // reduce to the form y^3 + py + q = 0
    // substitution: x = y - b/3a
    p = c/a - b*b/(3*a*a);
    q = 2*b*b*b/(27*a*a*a) - b*c/(3*a*a) + d/a;

    Q = p*p*p/27 + q*q/4;
    // Q > 0 - one real root and two complex conjugated roots
    // Q = 0 - one single real root and one double real root, or,
	//         if p = q = 0, then one triple real root
	// Q < 0 - three real roots

    if (Q >= 0) {
        //alpha   = mypowf(-q/2 + sqrtf(Q), 1.0/3);
        //beta    = mypowf(-q/2 - sqrtf(Q), 1.0/3);
        alpha   = cbrtf(-q/2 + sqrtf(Q));
        beta    = cbrtf(-q/2 - sqrtf(Q));
    } else {
        alpha   = cpowf(-q/2 + sqrtf(-Q) * I, 1.0/3);
        beta    = cpowf(-q/2 - sqrtf(-Q) * I, 1.0/3);
    }
    
//    cparts(alpha, &rp[0], &ip[0]);
//    cparts(beta, &rp[1], &ip[1]);

    roots[0] = alpha + beta - b/(3*a);
    roots[1] = -(alpha+beta)/2 - b/(3*a) + (alpha-beta)*sqrtf(3)/2 * I;
    roots[2] = -(alpha+beta)/2 - b/(3*a) - (alpha-beta)*sqrtf(3)/2 * I;
    
//    cparts(roots[0], &rp[2], &ip[2]);
//    cparts(roots[1], &rp[3], &ip[3]);
//    cparts(roots[2], &rp[4], &ip[4]);
    
    return 1;
}

/*
* solves equation ax^4 + bx^3 + cx^2 + dx + e = 0
* puts roots in the array
* returns 1 on success
*/
int solve4(float a, float b, float c, float d, float e, complexf roots[4])
{
    float p, q, r;
    float A, B, C, D, s;
    complexf cube[3];
    float a1, b1, c1, a2, b2, c2;
    int i;
    // float rp[3], ip[3]; // debug
    
    b /= a; c /= a; d /= a; e /= a;
    a = b; b = c; c = d; d = e;

    // reduce to the form y^4 + p*y^2 + q*y + r = 0
    p = b - 3*a*a/8;
    q = a*a*a/8 - a*b/2 + c;
    r = - 3*a*a*a*a/256 + a*a*b/16 - c*a/4 + d;

    // obtain cubic resolvent A*s^3 + B*s^2 + C*s + D = 0
    A = 2;
    B = -p;
    C = -2*r;
    D = r*p - q*q/4;
    if (!solve3(A, B, C, D, cube)) {
        return 0;
    }

//    for (i = 0; i < 3; i++) {
//        rp[i] = crealf(cube[i]);
//        ip[i] = cimagf(cube[i]);
//    }
    
    s = 0;
    if (crealf(cube[0]) > p/2) {
        s = crealf(cube[0]);
    } else if (crealf(cube[1]) > p/2) {
        s = crealf(cube[1]);
    } else if (crealf(cube[2]) > p/2) {
        s = crealf(cube[2]);
    }

    a1 = 1; b1 = -sqrtf(2*s-p); c1 = q/(2*sqrtf(2*s-p)) + s;
    a2 = 1; b2 = sqrtf(2*s-p); c2 = -q/(2*sqrtf(2*s-p)) + s;

    if (!solve2(a1, b1, c1, &roots[0]) || !solve2(a2, b2, c2, &roots[2])) {
        return 0;
    }
    for (i = 0; i < 4; i++) {
        roots[i] -= a/4;
    }
    return 1;
}

/*
* Calculates the distance from the point x to the ellipse with axis a, b
* *pd - output distance
* *pl - output auxilliary parameter lambda
* Returns 1 on success
*/
int dist_to_ellipse(float a, float b, float x[2], float *pd, float *pl)
{
    complexf roots[4];
    float d; // distance
    float l; // lambda
    int i;
    
    if (x[0] == 0 && x[1] == 0)
        return fminf(a, b);

    // Equation from necessary condition of conditional extrema
    // Looking for point on ellipse (Xe, Ye):
    // (x-Xe)^2 + (y-Ye)^2 -> min, when Xe^2/a^2 + Ye^2/b^2 = 1   
    if (!solve4(1.0f,
                2.0f * (a*a + b*b),
                a*a*a*a + b*b*b*b + 4*a*a*b*b - a*a*x[0]*x[0] - b*b*x[1]*x[1],
                2.0f * a*a*b*b * (a*a + b*b - x[0]*x[0] - x[1]*x[1]),
                a*a*b*b * (a*a*b*b - b*b*x[0]*x[0] - a*a*x[1]*x[1]),
                roots)) {
        return 0;
    }

    // Choose lambda that gives minimal distance
    *pd = FLT_MAX;
    for (i = 0; i < 4; i++) {
        if (cimagf(roots[i]) == 0) {
            l = crealf(roots[i]);
            
#define DIST_TO_ELLIPSE_EQ_PREC 1e-3f
            if (fabsf(a*a + l > DIST_TO_ELLIPSE_EQ_PREC) && 
                fabsf(b*b + l) > DIST_TO_ELLIPSE_EQ_PREC) 
            {
                d = sqrtf(powf(x[0]*l/(a*a + l), 2) + powf(x[1]*l/(b*b + l), 2));
                if (d < *pd) {
                    *pd = d;
                    *pl = l;
                }
            }
        }
    }
    return 1;
}

#define MAX_ITER_COUNT              50
#define MAX_ERROR_INCREASE_COUNT    4
#define PARAMS_COUNT                5
#define REL_PREC                    1e-5f

void global_to_canonical(float X[2], float alpha, float Xc[2], float x[2])
{
    float ca = cosf(alpha);
    float sa = sinf(alpha);
    x[0] = ca * (X[0] - Xc[1]) + sa * (X[1] - Xc[1]);
    x[1] = -sa * (X[0] - Xc[1]) + ca * (X[1] - Xc[1]);
}

void fill_jacobian( float *Ji, float x[2], float d, float l, float a, float b, float alpha)
{
    float dg_dl, dg_dx, dg_dy, dg_da, dg_db;
    float dl_dx, dl_dy, dl_da, dl_db;
    float dx_dXc, dx_dYc, dx_dalpha, dx_da, dx_db;
    float dy_dXc, dy_dYc, dy_dalpha, dy_da, dy_db;
    float dl_dXc, dl_dYc, dl_dalpha;
    float da2_dXc, da2_dYc, da2_da, da2_db, da2_dalpha;
    float db2_dXc, db2_dYc, db2_da, db2_db, db2_dalpha;
    
    if (l == 0 && d == 0) {
        memset(Ji, 0, PARAMS_COUNT * sizeof(float));
        return;
    }
    
    // auxiliary variables
    dg_dl = 4*l*l*l + 6*l*l * (a*a + b*a) + 
        2*l * (a*a*a*a + b*b*b*b + 4*a*a*b*b - a*a*x[0]*x[0] - b*b*x[1]*x[1]) + 
        2*a*a*b*a * (a*a + b*a - x[0]*x[0] - x[1]*x[1]);

    dg_dx = -2*x[0]*a*a * powf(b*b + l, 2);

    dg_dy = -2*x[1]*b*b * powf(a*b + l, 2);

    dg_da = l*l*l * 4*a + l*l * 2*a * (2*a*a + 4*b*b - x[0]*x[0]) + 
        l * 4*a*b*b * (2*a*a + b*b - x[0]*x[0] - x[1]*x[1]) + 
        2*a*b*b * (2*a*a*b*b - b*b*x[0]*x[0] - 2*a*a*x[1]*x[1]);

    dg_db = l*l*l * 4*b + l*l * 2*b * (2*b*b + 4*a*a - x[1]*x[1]) + 
        l * 4*a*a*b * (a*a + 2*b*b - x[0]*x[0] - x[1]*x[1]) + 
        2*a*a*b * (2*a*a*b*b - 2*b*b*x[0]*x[0] - a*a*x[1]*x[1]);
    
    
    dl_dx = - dg_dx / dg_dl;
    dl_dy = - dg_dy / dg_dl;
    dl_da = - dg_da / dg_dl;
    dl_db = - dg_db / dg_dl;
    
    dx_dXc = -cosf(alpha);
    dx_dYc = -sinf(alpha);
    dx_dalpha = x[1];
    dx_da = 0;
    dx_db = 0;
    
    dy_dXc = sinf(alpha);
    dy_dYc = -cosf(alpha);
    dy_dalpha = -x[0];
    dy_da = 0;
    dy_db = 0;
    
    dl_dXc = dl_dx * dx_dXc + dl_dy * dy_dXc;
    dl_dYc = dl_dx * dx_dYc + dl_dy * dy_dYc;
    dl_dalpha = dl_dx * dx_dalpha + dl_dy * dy_dalpha;
    
    da2_dXc = 0;
    da2_dYc = 0;
    da2_da = 2*a;
    da2_db = 0;
    da2_dalpha = 0;
    
    db2_dXc = 0;
    db2_dYc = 0;
    db2_da = 0;
    db2_db = 2*b;
    db2_dalpha = 0;
    
    // dd_dXc
    Ji[0] = (l/d) * (
        x[0] / powf(a*a + l, 3) * (a*a*x[0] * dl_dXc + l*(a*a + l) * dx_dXc - l*x[0] * da2_dXc)
        +
        x[1] / powf(b*b + l, 3) * (b*b*x[1] * dl_dXc + l*(b*b + l) * dy_dXc - l*x[1] * db2_dXc)
    );
    
    // dd_dYc
    Ji[1] = (l/d) * (
        x[0] / powf(a*a + l, 3) * (a*a*x[0] * dl_dYc + l*(a*a + l) * dx_dYc - l*x[0] * da2_dYc)
        +
        x[1] / powf(b*b + l, 3) * (b*b*x[1] * dl_dYc + l*(b*b + l) * dy_dYc - l*x[1] * db2_dYc)
    );
    
    // dd_da
    Ji[2] = (l/d) * (
        x[0] / powf(a*a + l, 3) * (a*a*x[0] * dl_da + l*(a*a + l) * dx_da - l*x[0] * da2_da)
        +
        x[1] / powf(b*b + l, 3) * (b*b*x[1] * dl_da + l*(b*b + l) * dy_da - l*x[1] * db2_da)
    );
    
    // dd_db
    Ji[3] = (l/d) * (
        x[0] / powf(a*a + l, 3) * (a*a*x[0] * dl_db + l*(a*a + l) * dx_db - l*x[0] * da2_db)
        +
        x[1] / powf(b*b + l, 3) * (b*b*x[1] * dl_db + l*(b*b + l) * dy_db - l*x[1] * db2_db)
    );
    
    // dd_alpha
    Ji[4] = (l/d) * (
        x[0] / powf(a*a + l, 3) * (a*a*x[0] * dl_dalpha + l*(a*a + l) * dx_dalpha - l*x[0] * da2_dalpha)
        +
        x[1] / powf(b*b + l, 3) * (b*b*x[1] * dl_dalpha + l*(b*b + l) * dy_dalpha - l*x[1] * db2_dalpha)
    );
}

/*
* LU decomposition of matrix A
* Matrices L and U are return in the respective parts of the matrix A 
* indx - row perturbations array
* Return 1 on success, 0 if the matrix is singular
*/
int ludcmp(float **a, int n, int *indx) 
{
    int i, imax, j, k;
    float big, dum, sum, temp;
    float vv[MAX_DCMP_N];

    // calculate the biggest elements in every row
    for (i = 0; i < n; i++) {
        big = 0.0f;
        for (j = 0; j < n; j++) { 
            if ((temp = fabsf(a[i][j])) > big)
                big = temp;
        }
        if (big == 0.0f)
            return 0; // Singular matrix in routine ludcmp
        vv[i] = (1.0f / big);
    }

    for (j = 0; j < n; j++) {
        // calculate elements of U
        for (i = 0; i < j; i++) {
            sum = a[i][j];
            for (k = 0; k < i; k++) 
                sum -= (a[i][k] * a[k][j]);
            a[i][j] = sum;
        }
        // calculate elements of L
        big = 0.0f;
        for (i = j; i < n; i++) {
            sum = a[i][j];
            for (k = 0; k < j; k++)
                sum -= a[i][k] * a[k][j];
            a[i][j] = sum;

            // find the row with the biggest jth element
            if ((dum = vv[i] * fabsf(sum)) >= big) {
                big = dum;
                imax = i;
            }
        }
        // figure out if we need to interchange rows
        if (j != imax) {
            for (k = 0; k < n; k++) {
                dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dum;
            }
            vv[imax] = vv[j];
        }
        // remember where the jth row goes
        indx[j] = imax;

        // carry out division by the revealed biggest jj-th element
        dum = 1.0f / a[j][j];
        for (i = j + 1; i < n; i++)
            a[i][j] *= dum;
    }
    //free(vv);
    return 1;
}

/*
* Solves the set of linear equations Ax = b
* a     - matrix LU-decomposition returned by ludcmp
* indx  - vector of row perturbations returned by ludcmp
* b     - right side vector and solution is returned in it 
*/
void lubksb(float **a, int n, float *b, int *indx)
{
    int i, ii = -1, ip, j;
    float sum;

    for (i = 0; i < n; i++) {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if (ii >= 0) {
            for (j = ii; ii < i; ii++)
                sum -= a[i][j] * b[j];
        }
        else if (sum)
            ii = i;
        b[i] = sum;
    }
    for (i = n - 1; i >= 0; i--) {
        sum = b[i];
        for (j = i + 1; j < n; j++)
            sum -= a[i][j] * b[j];
        b[i] = sum / a[i][i];
    }
}

/*
* Solves system of linear algebraic equations Ax = b
* Returns the solution in matrix b
* Returns 1 on success
*/
int linsolve(arm_matrix_instance_f32 *pA, arm_matrix_instance_f32 *pb)
{
    float *a[MAX_DCMP_N], *pdata;
    int i, n, indx[MAX_DCMP_N];
    
    n = pA->numRows;
    pdata = pA->pData;
    
    for (i = 0; i < pA->numRows; i++) {
        a[i] = &(pdata[i*n]);
    }
    
    if (!ludcmp(a, n, indx))
        return 0;
    
    lubksb(a, n, pb->pData, indx);
    return 1;
}

// check potentially zero parameter
uint8_t check_pot_zero_param(float value, float delta)
{
    if (fabsf(value) < REL_PREC) {  // == 0
        return (abs(delta) < REL_PREC) ? 1 : 0;
    } else {
        return (abs(delta/value) < REL_PREC) ? 1 : 0;
    }
}

float norm(float *v, int n) 
{
    float result;
    arm_dot_prod_f32(v, v, n, &result);
    return sqrtf(result);
}

int ellipse_fitting(float **points, int N, 
                    float init_Xc[2], float init_a, float init_b, float init_alpha,
                    float Xc[2], float *pa, float *pb, float *palpha)
{
    int i;
    int iter_count = 1;
    float a, b, alpha;
    float x[2], d, l;
    
    float *J_data, *Jt_data, *e_data;
    float Js_data[PARAMS_COUNT * PARAMS_COUNT];
    float es_data[PARAMS_COUNT];
    arm_matrix_instance_f32 J, Jt, e, Js, es, da;
    
    uint8_t Xc_ok[2], alpha_ok;
    uint8_t rel_prec_achieved;
    
    Xc[0]   = init_Xc[0]; 
    Xc[1]   = init_Xc[1];
    a       = init_a;
    b       = init_b;
    alpha   = init_alpha;
    
    // alloc memory 
    J_data  = (float *)malloc(N * PARAMS_COUNT * sizeof(float));
    Jt_data = (float *)malloc(N * PARAMS_COUNT * sizeof(float));
    e_data  = (float *)malloc(N * sizeof(float));
    
    while (iter_count <= MAX_ITER_COUNT) {
        for (i = 0; i < N; i++) {
            global_to_canonical(points[i], alpha, Xc, x);
            if (!dist_to_ellipse(a, b, x, &d, &l)) {
                return 0;
            }
            e_data[i] = -d;
            
            
            if (fabsf(1 - fabsf(a/b)) < REL_PREC) {
                fill_jacobian(J_data + i*(PARAMS_COUNT-1), x, d, l, a, b, alpha);
            } else {
                fill_jacobian(J_data + i*PARAMS_COUNT, x, d, l, a, b, alpha);
            }
        }
        
        // solve linear system J * da = e
        // to find ellipse parameters' changes
        if (fabsf(1 - fabsf(a/b)) < REL_PREC) {
            arm_mat_init_f32(&J, N, PARAMS_COUNT-1, J_data);
            arm_mat_init_f32(&Jt, PARAMS_COUNT-1, N, Jt_data);
            arm_mat_init_f32(&e, N, 1, e_data);
            arm_mat_init_f32(&Js, PARAMS_COUNT-1, PARAMS_COUNT-1, Js_data);
            arm_mat_init_f32(&es, PARAMS_COUNT-1, 1, es_data);

            arm_mat_trans_f32(&J, &Jt);
            arm_mat_mult_f32(&Jt, &J, &Js);
            arm_mat_mult_f32(&Jt, &e, &es);

            if (!linsolve(&Js, &es)) {
                return 0;
            }
        } else {
            arm_mat_init_f32(&J, N, PARAMS_COUNT, J_data);
            arm_mat_init_f32(&Jt, PARAMS_COUNT, N, Jt_data);
            arm_mat_init_f32(&e, N, 1, e_data);
            arm_mat_init_f32(&Js, PARAMS_COUNT, PARAMS_COUNT, Js_data);
            arm_mat_init_f32(&es, PARAMS_COUNT, 1, es_data);

            arm_mat_trans_f32(&J, &Jt);
            arm_mat_mult_f32(&Jt, &J, &Js);
            arm_mat_mult_f32(&Jt, &e, &es);

            if (!linsolve(&Js, &es)) {
                return 0;
            }
            alpha += es_data[4];
        }
        // solution is in the vector ls
        Xc[0]   += es_data[0];
        Xc[1]   += es_data[1];
        a       += es_data[2];
        b       += es_data[3];
        
        memset(Xc_ok, 0, 2 * sizeof(uint8_t));
        for (i = 0; i < 2; i++) {
            Xc_ok[i] = check_pot_zero_param(Xc[i], es_data[i]);
        }
        
        alpha_ok = 0;
        if (da.numRows < 5) 
            alpha_ok = 1;
        else if (check_pot_zero_param(alpha, es_data[4]))
            alpha_ok = 1;

        rel_prec_achieved = Xc_ok[0] && Xc_ok[1] && alpha_ok && 
                            fabsf(es_data[2]/a) < REL_PREC && fabs(es_data[3]/b) < REL_PREC;
        
        if (rel_prec_achieved || norm(es_data, da.numRows) < powf(REL_PREC, 1.5)) {
            *pa = a;
            *pb = b;
            *palpha = alpha;
            return 1;
        }
        
        iter_count++;
    }
    return 0;
}
