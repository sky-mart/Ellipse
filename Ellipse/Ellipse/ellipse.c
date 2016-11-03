#include "ellipse.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <mkl.h>
#include <time.h>
#include <math.h>

static const int SPACE_SIZE		= 2; // space dimenstion
static const int PARAMS_COUNT	= 5; // ellipse parameters number

void printSystem(myreal *A, myreal *b, int m, int n)
{
	int i, j;
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			printf("%.4lf\t", A[i*n + j]); 
		}
		printf("\t%.4lf\n", b[i]);
	}
}

void checkJandE(myreal Ji[2][5], myreal e[2])
{
	int i, j;
	for (i = 0; i < 2; i++) {
		for (j = 0; j < 5; j++) {
			if (abs(Ji[i][j]) > 1e+3) {
				printf("J Alarm!\n");
			}
		}
		if (e[i] > 1) {
			printf("E Alarm!\n");
		}
	}
}

void ellipseFitting(myreal **X, int size, myreal ellipseParams[5])
{
	int iter = 0;				// iteration counter
	myreal *J;					// jacobian
	myreal *e;					// penalty vector
	myreal x[5];				// solution of linear eqauations system
	myreal ea, eb, s;			// for loop break conditions
	int i;		
	
	// initial parameters for the algorithm are taken from the mass center
	massCenter(X, size, ellipseParams);
	ellipseParams[3] = ellipseParams[2];
	ellipseParams[4] = 0;
	/*ellipseParams[0] = -7.000022;
	ellipseParams[1] = 21.000046;
	ellipseParams[2] = 129.819432;
	ellipseParams[3] = 90.102916;
	ellipseParams[4] = 0.394887;*/

	// memory allocation
	e = (myreal *)malloc(size * SPACE_SIZE * sizeof(myreal));
	J = (myreal *)malloc(size * SPACE_SIZE * PARAMS_COUNT * sizeof(myreal));

	do {	
		myreal Ji[2][5];
		myreal ei[2];
		myreal start, finish, average = 0;
		//start = clock();
		for (i = 0; i < size; i++) {
			// get jacobian part corresponding to ith point
			//printf("%i ", i);
			jacobianAndPenalty(X[i], ellipseParams, Ji, ei);
			//checkJandE(Ji, ei);
		
			// copy the part to the whole
			memcpy(J + i*PARAMS_COUNT*SPACE_SIZE, Ji, PARAMS_COUNT*SPACE_SIZE*sizeof(myreal));
			memcpy(e + i*SPACE_SIZE, ei, SPACE_SIZE*sizeof(myreal));
		}
		//finish = clock();
		//average = (finish - start) / size / CLOCKS_PER_SEC;
		//printf("Average J calculation: %.12lf\n", average);

		//printSystem(J, e, size*2, 5);
		if (fabs(ellipseParams[2] - ellipseParams[3]) < 1e-5) { // a == b
			linsolve(J, e, x, size*SPACE_SIZE, PARAMS_COUNT - 1);
			x[4] = 0;
		} else {
			linsolve(J, e, x, size*SPACE_SIZE, PARAMS_COUNT);
		}

		//printf("popravki:\n");
		for (i = 0; i < 5; i++) {
			//printf("%i: %lf - %lf\n", i, ellipseParams[i], x[i]);
			ellipseParams[i] += x[i];
		}

		ea = fabs(x[2] / ellipseParams[2]);
		eb = fabs(x[3] / ellipseParams[3]);
		s = sqrt(x[0]*x[0] + x[1]*x[1]);
		printf("ea=%lf, eb=%lf, s=%lf\n", ea, eb, s);

		iter++;
	} while((ea >= 1e-3) || (eb >= 1e-3) ||
		(2*s/(ellipseParams[2]+ellipseParams[3]) >= 1e-3));

	//printf("%i iterations\n", iter);
	free(e);
	free(J);
}

void linsolve(myreal *J, myreal *e, myreal *x, int nRows, int nParams)
{
	myreal *Jf;		// matrix of exact size we need
	myreal *Jft;	// the same but transposed
	myreal *Js;		// square matrix
	myreal *Jst;	// the same but transposed
	int i, j;

	// parameters for gesv
	int nrhs;
	int ipiv[5];
	int lda;
	int ldb;
	int info;

	// adjust size
	if (nParams == PARAMS_COUNT) {
		Jf = J;
	} else {
		Jf = (myreal *)malloc(nRows * nParams * sizeof(myreal));
		for (i = 0; i < nRows; i++) {
			memcpy(Jf + i*nParams, J + i*PARAMS_COUNT, nParams*sizeof(myreal));
		}
	}

	// squarify our matrix to get rid of overdetermination
	Js = (myreal *)malloc(nParams * nParams * sizeof(myreal));
	cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, nParams, nParams, nRows,
				1, Jf, nParams, Jf, nParams, 0, Js, nParams);

	// transpose for gemv
	Jft = (myreal *)malloc(nRows * nParams * sizeof(myreal));
	for (i = 0; i < nParams; i++) {
		for (j = 0; j < nRows; j++) {
			Jft[i*nRows + j] = Jf[j*nParams + i];
		}
		x[i] = 0;
	}

	// decrease the right side column size in appropriate way
	// put it in x
	cblas_dgemv(CblasRowMajor, CblasNoTrans, nParams, nRows, 
		1, Jft, nRows, e, 1, 0, x, 1);

	// transpose for gesv
	Jst = (myreal *)malloc(nParams * nParams * sizeof(myreal));
	for (i = 0; i < nParams; i++) {
		for (j = 0; j < nParams; j++) {
			*(Jst + i*nParams + j) = *(Js + j*nParams + i);
		}
	}

	// solve linear equations system
	nrhs = 1; 
	lda = nParams; 
	ldb = nParams; 
	dgesv(&nParams, &nrhs, Jst, &lda, ipiv, x, &ldb, &info);

	if (nParams != PARAMS_COUNT) {
		free(Jf);
	}
	free(Jft);
	free(Js);
	free(Jst);
}

// Ji[2][5];	part of the jacobian for the point
// ei[2];		penalty for this point
void jacobianAndPenalty(myreal Xi[2], myreal ellipseParams[5], myreal Ji[2][5], myreal ei[2])
{
	myreal Xc, Yc, a, b, alpha;				// ellipse parameters
	myreal C, S, R[2][2], invR[2][2];		// rotation matrix and its composite parts
	myreal Q[2][2], invQ[2][2], M[2][5];	// auxiliary matrices
	myreal inter[2][5];

	myreal x[2];	// point coordinates in the ellipse coordinate system
	myreal xs[2];	// the closest point coordinates in the ellipse coordinate system
	myreal Xs[2];	// the closest point coordinates in the global coordinate system

	Xc		= ellipseParams[0];
	Yc		= ellipseParams[1];
	a		= ellipseParams[2];
	b		= ellipseParams[3];
	alpha	= ellipseParams[4];

	C = cos(alpha);
	S = sin(alpha);
	R[0][0] = C;  R[0][1] = S;
	R[1][0] = -S; R[1][1] = C;

	conversion(Xi, R, ellipseParams, x);		// transition to the ellipse coordinate system
	distToEllipse(a, b, x[0], x[1], xs);		// the closest point search
	invConversion(xs, R, ellipseParams, Xs);	// transition to the global coordinate system

	// form auxiliary matrices
	Q[0][0] = b*b * xs[0];
	Q[0][1] = a*a * xs[1];
	Q[1][0] = (a*a - b*b)*xs[1] + b*b*x[1];
	Q[1][1] = (a*a - b*b)*xs[0] - a*a*x[0];

	B(a, b, C, S, xs, x, M);
	simpleInv(Q, invQ);
	simpleInv(R, invR);

	// Ji = inv(R) * inv(Q) * M; 

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, SPACE_SIZE, PARAMS_COUNT, SPACE_SIZE,
		1, invQ, SPACE_SIZE, M, PARAMS_COUNT, 0, inter, PARAMS_COUNT);

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, SPACE_SIZE, PARAMS_COUNT, SPACE_SIZE,
		1, invR, SPACE_SIZE, inter, PARAMS_COUNT, 0, Ji, PARAMS_COUNT);

	// penalty ei = X + al*Xs
	cblas_dcopy(SPACE_SIZE, Xi, 1, ei, 1);
	cblas_daxpy(SPACE_SIZE, -1, Xs, 1, ei, 1); 
}

void simpleInv(myreal A[2][2], myreal invA[2][2])
{
	myreal det = A[0][0]*A[1][1] - A[0][1]*A[1][0];
	invA[0][0] = A[1][1] / det;
	invA[0][1] = -A[0][1] / det;
	invA[1][0] = -A[1][0] / det;
	invA[1][1] = A[0][0] / det;
}

void B(myreal a, myreal b, myreal C, myreal S, myreal R[2], myreal Ri[2], myreal M[2][5])
{
	M[0][0] = b*b*R[0]*C - a*a*R[1]*S;
    M[1][0] = b*b*(Ri[1]-R[1])*C + a*a*(Ri[0]-R[0])*S;

	M[0][1] = b*b*R[0]*S + a*a*R[1]*C;
    M[1][1] = b*b*(Ri[1]-R[1])*S - a*a*(Ri[0]-R[0])*C;

	M[0][2] = a*(b*b - R[1]*R[1]);
    M[1][2] = 2*a*R[1]*(Ri[0] - R[0]);

	M[0][3] = b*(a*a - R[0]*R[0]);
    M[1][3] = -2*b*R[0]*(Ri[1] - R[1]);

	M[0][4] = (a*a - b*b)*R[0]*R[1];
    M[1][4] = (a*a - b*b)*(R[0]*R[0] - R[1]*R[1] - R[0]*Ri[0] + R[1]*Ri[1]);
}

void massCenter(myreal **X, int n, myreal params[3])
{
	int i;
	myreal Xc[2], R;
	Xc[0] = 0; Xc[1] = 0;
	for (i = 0; i < n; i++) { 
		cblas_daxpy(2, 1, X[i], 1, Xc, 1);
	}
	cblas_dscal(2, 1.0 / n, Xc, 1);
    
    R = 0;
    for (i = 0; i < n; i++) { 
        R += pow(Xc[0] - X[i][0], 2) + pow(Xc[1] - X[i][1], 2);
	}
    R = sqrt(R/n);
	
	params[0] = Xc[0];
	params[1] = Xc[1];
	params[2] = R;
}

void conversion(myreal Rglob[2], myreal T[2][2], myreal C[2], myreal Rel[2])
{
	myreal tmp[2];
	cblas_dcopy(SPACE_SIZE, Rglob, 1, tmp, 1);
	cblas_daxpy(SPACE_SIZE, -1, C, 1, tmp, 1);
	cblas_dgemv(CblasRowMajor, CblasNoTrans, SPACE_SIZE, SPACE_SIZE, 1, T, SPACE_SIZE, 
		tmp, 1, 0, Rel, 1);
}

void invConversion(myreal Rel[2], myreal T[2][2], myreal C[2], myreal Rglob[2])
{
	//Rglob = inv(R)*Rel + C;
	myreal tmp[2];
	myreal invT[2][2];
	simpleInv(T, invT);
	cblas_dcopy(SPACE_SIZE, Rel, 1, tmp, 1);
	cblas_dgemv(CblasRowMajor, CblasNoTrans, SPACE_SIZE, SPACE_SIZE, 1, invT, SPACE_SIZE, 
		tmp, 1, 0, Rglob, 1);
	cblas_daxpy(SPACE_SIZE, 1, C, 1, Rglob, 1);
}

void distToEllipse(myreal a, myreal b, myreal X0, myreal Y0, myreal d[2])
{
	myreal a0, a1, a2, a3, a4; // equation coefficients
	mycomplex roots[4];
	myreal minX, minY, minVal;
	int i;

	if (X0 == 0 && Y0 == 0) {
		
	} else {
		myreal m = a > b ? a : b;
		a /= m;
		b /= m;
		X0 /= m;
		Y0 /= m;

		// solve equation a0*l^4 + a1*l^3 + a2*l^2 + a3*l + a4 = 0	
		a0 = - 1 / (a*a * b*b);
		a1 = - 2/(a*a) - 2/(b*b);
		a2 = X0*X0/(b*b) + Y0*Y0/(a*a) - b*b/(a*a) - a*a/(b*b) - 4;
		a3 = 2*X0*X0 + 2*Y0*Y0 - 2*b*b - 2*a*a;
		a4 = X0*X0*b*b + Y0*Y0*a*a - a*a*b*b;

		solve4(a0, a1, a2, a3, a4, roots,a,b);

		minX = -3.14e+28;
		minY = -3.14e+28;
		minVal = pow(minX - X0, 2) + pow(minY - Y0, 2);

		for (i = 0; i < 4; i++) {
			myreal absIm = roots[i].i >= 0 ? roots[i].i : -roots[i].i;
			if (absIm < 1e-15) {//((absIm < 1e-10) || (abs(roots[i].i/roots[i].r) < 1e-10)) {
				myreal root = roots[i].r;
				myreal tmpX, tmpY;
				if (((2 + 2*root/a/a > 0) && (2 + 2*root/b/b >= 0)) || 
						((2 + 2*root/a/a >= 0) && (2 + 2*root/b/b > 0))) {
					if (root == -a*a) {
						tmpY = Y0 / (1 + root/b/b);
						tmpX = sign(tmpY) * (b/a) * sqrt(a*a - tmpY*tmpY);
					} else if (root == -b*b) {
						tmpX = X0 / (1 + root/a/a);
						tmpY = sign(tmpX) * (a/b) * sqrt(b*a - tmpX*tmpX);
					} else {
						tmpX = X0 / (1 + root/a/a);
						tmpY = Y0 / (1 + root/b/b);
					}

					if ((tmpX * X0 >= 0) && (tmpY * Y0 >= 0)) {
						myreal val = pow(tmpX - X0, 2) + pow(tmpY - Y0, 2);
						if (val < minVal) {
							minX = tmpX;
							minY = tmpY;
							minVal = val;
						}
					}
				}
			}
		}
		d[0] = minX * m;
		d[1] = minY * m;
	}
} 

void solve4(myreal a, myreal b, myreal c, myreal d, myreal e, mycomplex roots[4])
// Solve the 4th degree equation x^4 + a*x^3 + b*x^2 + c*x + d = 0
{
	int j,k,i;
	myreal p, q, r;					// coefficients in the 3rd degree equation
	myreal A, B, C, D;				// coefficients in cubic resolvent
	mycomplex cubeRoots[3];			// cubic resolvent roots
	mycomplex iRoots[4];
	mycomplex s;					// auxuliary expression
	mycomplex sqcf1[3], sqcf2[3];
	mycomplex sqr1[2], sqr2[2];		// square equations roots
	myreal cErr;
	myreal tErr;
	mycomplex tmp;
	b = b/a; c = c/a; d = d/a; e = e/a;
	a = b; b = c; c = d; d = e;

	// here is the form without a 3rd degree
	// y^4 + p*y^2 + q*y + r = 0
	p = b - 3*a*a/8;
	q = a*a*a/8 - a*b/2 + c;
	r = - 3*a*a*a*a/256 + a*a*b/16 - c*a/4 + d;

	// Obtain cubic resolvent
	// A*s^3 + B*s^2 + C*s + D = 0
	A = 2;
	B = -p;
	C = -2*r;
	D = r*p - q*q/4;

	solve3(A, B, C, D, cubeRoots);
	s = makeComplex(0, 0);
	if (cubeRoots[0].r > p/2) { 
		s = cubeRoots[0];
	} else if (cubeRoots[1].r > p/2) {
		s = cubeRoots[1];
	} else if (cubeRoots[2].r > p/2) {
		s = cubeRoots[2];
	}

	sqcf1[0] = makeComplex(1, 0);
	sqcf1[1] = zmulr(complexNthroot(zsubr(zmulr(s, 2), p), 2), -1);
	sqcf1[2] = zaddz(zdivz(makeComplex(q, 0), zmulr(complexNthroot(zsubr(zmulr(s, 2), p), 2), 2)),s);
	solve2(sqcf1[0], sqcf1[1], sqcf1[2], sqr1);

	sqcf2[0] = makeComplex(1, 0);
	sqcf2[1] = complexNthroot(zsubr(zmulr(s, 2), p), 2);
	sqcf2[2] = zaddz(zdivz(makeComplex(-q, 0), zmulr(complexNthroot(zsubr(zmulr(s, 2), p), 2), 2)),s);
	solve2(sqcf2[0], sqcf2[1], sqcf2[2], sqr2);
		
	roots[0] = zsubr(sqr1[0], a/4);
	roots[1] = zsubr(sqr1[1], a / 4);
	roots[2] = zsubr(sqr2[0], a / 4);
	roots[3] = zsubr(sqr2[1], a / 4);
}

void solve3(myreal a, myreal b, myreal c, myreal d, mycomplex roots[3])
// Solve cubic equation ax^3 + bx^2 + cx + d = 0
{
	myreal p, q;		// coefficients in the equation without a 2nd degree
	myreal Q;			// cubic equation determinant
	mycomplex al, be;	// auxiliary expressions
	mycomplex y[3];		// reduced equation roots 
	myreal tmp;
	myreal prec = 1e-12;

	// get form  y^3 + py + q = 0
	// replacement: x = y - b/3a
	p = c/a - b*b/(3*a*a);
	q = 2*b*b*b/(27*a*a*a) - b*c/(3*a*a) + d/a;

	Q = pow(p/3, 3) + pow(q/2, 2);
	// Q > 0 - one real root and two complex conjugated roots
	// Q = 0 - one single real root and one double real root, or,
	//         if p = q = 0, then one triple real root
	// Q < 0 - three real roots

	if (Q >= 0) { // in this case everything is fine and we simply extract root from a real number
		tmp = -q/2 + sqrt(Q);
		//al.r = sign(tmp) * pow(abs(tmp), 1/3.0);//exp(log(abs(tmp)) / 3);
		al.r = cbrt(tmp); //cbrt(tmp, prec);
		al.i = 0;
		tmp = -q/2 - sqrt(Q);
		//be.r = sign(tmp) * pow(abs(tmp), 1/3.0);//exp(log(abs(tmp)) / 3);
		be.r = cbrt(tmp); //cbrt(tmp, prec);
		be.i = 0;
	} else { // here we have to compose a complex number: Re = -q/2, Im = +-sqrt(Q), 
		// then, choose the first root variant out of three
		al = complexNthroot(makeComplex(-q/2, sqrt(-Q)), 3);
		be = complexNthroot(makeComplex(-q/2, -sqrt(-Q)), 3);
	}
	
	y[0] = zaddz(al, be);
	y[1] = zaddz(zdivr(zaddz(al, be), -2), zmulz(zsubz(al, be), makeComplex(0, sqrt(3)/2))); 
	y[2] = zsubz(zdivr(zaddz(al, be), -2), zmulz(zsubz(al, be), makeComplex(0, sqrt(3)/2))); 

	roots[0] = zsubr(y[0], b/(3*a));
	roots[1] = zsubr(y[1], b/(3*a));
	roots[2] = zsubr(y[2], b/(3*a));
}

void solve2(mycomplex a, mycomplex b, mycomplex c, mycomplex roots[2])
{
	mycomplex D;
	
	if (zabs(a) == 0) {
		if (zabs(b) == 0) return;
		roots[0] = roots[1] = zmulr(zdivz(c, b), -1);
	}

	D = zsubz(zmulz(b, b), zmulr(zmulz(a, c), 4));
	if (zabs(D) == 0) {
		roots[0] = roots[1] = zdivz(zmulr(b, -1), zmulr(a, 2));
	}  else {
		roots[0] = zdivz(zsubz(zmulr(b, -1), complexNthroot(D, 2)), zmulr(a, 2));
		roots[1] = zdivz(zaddz(zmulr(b, -1), complexNthroot(D, 2)), zmulr(a, 2));
	}
}

int nDigits(myreal number)
{
	int intNumber = (int)ceil(number);
	int count = 0;
	while (intNumber > 0) {
		intNumber /= 10;
		count++;
	}
	return count;
}

//myreal cbrt(myreal a, myreal e)
//{
//	myreal x[2];
//	int n;
//
//	if (a < 0) {
//		return - cbrt(-a, e);
//	}
//
//	n = nDigits(a);
//	x[1] = pow(10, n / 3);
//	do {
//		x[0] = x[1];
//		x[1] = 2*x[0]/3 + a/3/x[0]/x[0];
//	} while ((x[1]*x[1]*x[1] - a) >= e);
//	return x[1];
//}

