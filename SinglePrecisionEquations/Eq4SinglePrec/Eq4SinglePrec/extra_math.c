//
//  extra_math.c
//  Eq4SinglePrec
//
//  Created by Влад Агиевич on 28.02.17.
//  Copyright © 2017 Mocsmart. All rights reserved.
//

#include "extra_math.h"
#include <stdio.h>
#include <stdlib.h>

void cparts(complexr z, real *pr, real *pi)
{
    *pr = crealr(z);
    *pi = crealr(z);
}

void cprint(complexr z)
{
    printf("%e + %ei", crealr(z), cimagr(z));
}

//real mypowr(real b, real e)
//{
//    return (b >= 0) ? powr(b, e) : -powr(-b, e);
//}

void dotprodf(float *a, float *b, uint n, float *pResult)
{
    uint i;
    *pResult = 0;
    for (i = 0; i < n; i++) {
        *pResult += a[i] * b[i];
    }
}

void dotprod(double *a, double *b, uint n, double *pResult)
{
    uint i;
    *pResult = 0;
    for (i = 0; i < n; i++) {
        *pResult += a[i] * b[i];
    }
}

void vec_sub(real *a, real *b, real *d, uint n)
{
    uint i;
    for (i = 0; i < n; i++) {
        d[i] = a[i] - b[i];
    }
}

void vec_maxf(float *v, uint n, float *pResult, uint *pIndex)
{
    uint i, bigIndex = 0;
    float big = v[0];
    for (i = 1; i < n; i++) {
        if (v[i] > big) {
            bigIndex = i;
            big = v[i];
        }
    }
    *pResult = big;
    *pIndex = bigIndex;
}

void vec_max(double *v, uint n, double *pResult, uint *pIndex)
{
    uint i, bigIndex = 0;
    double big = v[0];
    for (i = 1; i < n; i++) {
        if (v[i] > big) {
            bigIndex = i;
            big = v[i];
        }
    }
    *pResult = big;
    *pIndex = bigIndex;
}

float normf(float *v, uint n)
{
    float result;
    dotprodf(v, v, n, &result);
    return sqrtf(result);
}

double norm(double *v, uint n)
{
    double result;
    dotprod(v, v, n, &result);
    return sqrt(result);
}

real dist(real *u, real *v, int n)
{
    real val;
    real *d = (real *)malloc(n * sizeof(real));
    
    vsub(u, v, d, n);
    val = normr(d, n);
    free(d);
    return val;
}

void mat_init(matrix *pA, uint16_t nRows, uint16_t nCols, mreal *pData)
{
    pA->numRows = nRows;
    pA->numCols = nCols;
    pA->pData   = pData;
}

int mat_trans(matrix *pA, matrix *pAt)
{
    int i, j;
    for (i = 0; i < pA->numRows; i++) {
        for (j = 0; j < pA->numCols; j++) {
            pAt->pData[j * pAt->numCols + i] = pA->pData[i * pA->numCols + j];
        }
    }
    return 1;
}

int mat_sub(matrix *pA, matrix *pB, matrix *pC)
{
    int i, j;
    for (i = 0; i < pA->numRows; i++) {
        for (j = 0; j < pA->numCols; j++) {
            pC->pData[i * pC->numCols + j] = pA->pData[i * pA->numCols + j] -
            pB->pData[i *pB->numCols + j];
        }
    }
    return 1;
}

int mat_mult(matrix *pA, matrix *pB, matrix *pC)
{
    int i, j, k;
    mreal tmp;
    for (i = 0; i < pC->numRows; i++) {
        for (j = 0; j < pC->numCols; j++) {
            tmp = 0;
            for (k = 0; k < pA->numCols; k++) {
                tmp += pA->pData[i * pA->numCols + k] * pB->pData[k * pB->numCols + j];
            }
            pC->pData[i * pC->numCols + j] = tmp;
        }
    }
    return 1;
}

/*
 * LU decomposition of matrix A
 * Matrices L and U are return in the respective parts of the matrix A
 * indx - row perturbations array
 * Return 1 on success, 0 if the matrix is singular
 */
int ludcmp(mreal **a, int n, int *indx)
{
    int i, imax = 0, j, k;
    mreal big, dum, sum, temp;
    mreal vv[MAX_DCMP_N];
    
    // calculate the biggest elements in every row
    for (i = 0; i < n; i++) {
        big = 0.0;
        for (j = 0; j < n; j++) {
            if ((temp = fabsm(a[i][j])) > big)
                big = temp;
        }
        if (big == 0.0)
            return 0; // Singular matrix in routine ludcmp
        vv[i] = (1.0 / big);
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
        big = 0.0;
        for (i = j; i < n; i++) {
            sum = a[i][j];
            for (k = 0; k < j; k++)
                sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
            
            // find the row with the biggest jth element
            if ((dum = vv[i] * fabsm(sum)) >= big) {
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
        dum = 1.0 / a[j][j];
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
void lubksb(mreal **a, int n, mreal *b, int *indx)
{
    int i, ii = -1, ip, j;
    mreal sum;
    
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
int linsolve(matrix *pA, matrix *pb)
{
    mreal *a[MAX_DCMP_N], *pdata;
    int i, n, indx[MAX_DCMP_N];
    
    n = pA->numRows;
    pdata = pA->pData;
    
    for (i = 0; i < n; i++) {
        a[i] = &(pdata[i*n]);
    }
    
    if (!ludcmp(a, n, indx))
        return 0;
    
    lubksb(a, n, pb->pData, indx);
    return 1;
}
