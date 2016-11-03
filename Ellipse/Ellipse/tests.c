#include <stdio.h>

#include "tests.h"
#include "mycomplex.h"
#include "ellipse.h"
#include "mkl.h"
#include <stdlib.h>
#include <time.h>

int complexTest(int output)
{
	int res = 1;
	res = res && addTest(output);
	if (output) {
		if (res) {
			printf("complexTest: ALL TESTS PASSED\n");
		} else {
			printf("complexTest: SOME TESTS FAILED\n");
		}
	}
	return res;
}

int addTest(int output)
{
	int res = 0;
	mycomplex z1 = makeComplex(1, 4);
	mycomplex z2 = makeComplex(2, -8);
	mycomplex z = zaddz(z1,z2);
	res = (z.r == 3 && z.i == -4) ? 1 : 0;
	if (output) {
		if (res) {
			printf("addTest: PASSED\n");
		} else {
			printf("addTest: FAILED\n");
		}
	}
	return res;
}

int blasTest(int output)
{
	int res = 1;
	res = res && copyTest(output);
	res = res && axpyTest(output);
	res = res && gemvTest(output);
	res = res && gemmTest(output);
	if (output) {
		if (res) {
			printf("blasTest: ALL TESTS PASSED\n");
		} else {
			printf("blasTest: SOME TESTS FAILED\n");
		}
	}
	return res;
}

int copyTest(int output)
{
	int res = 1;
	int size = 4;
	myreal a[] = {1, 2, 3, 4};
	myreal b[] = {2, 4, 0, -1};
	int i;

	cblas_dcopy(size, a, 1, b, 1);
	
	for (i = 0; i < size; i++) {
		if (b[i] != a[i]) {
			res = 0;
			break;
		}
	}

	if (output) {
		if (res) {
			printf("copyTest: PASSED\n");
		} else {
			printf("copyTest: FAILED\n");
		}
	}
	return res;
}

int axpyTest(int output)
{
	int res = 1;
	int size = 4;
	myreal a[] = {1, 2, 3, 4};
	myreal b[] = {2, 4, 0, -1};
	myreal al = 1;
	myreal be = 0;

	cblas_daxpy(size, al, a, 1, b, 1);
	res == res && (b[0] == 3);
	res == res && (b[1] == 6);
	res == res && (b[2] == 3);
	res == res && (b[3] == 3);
	
	if (output) {
		if (res) {
			printf("axpyTest: PASSED\n");
		} else {
			printf("axpyTest: FAILED\n");
		}
	}
	return res;
}



int gemvTest(int output)
{
	int res = 1;
	myreal a[3][2] = {1, 2, 3, 4, 5, 6};
	myreal at[2][3] = {1, 3, 5, 2, 4, 6};
	myreal b[3] = {7, -1, 5};
	myreal c[2] = {0, 0};
	int m = 2, n = 3, k = 3;
	cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, 
   		1, at, n, b, 1, 0, c, 1);

	res = res && (c[0] == 29);
	res = res && (c[1] == 40);

	if (output) {
		if (res) {
			printf("gemvTest: PASSED\n");
		} else {
			printf("gemvTest: FAILED\n");
		}
	}
	return res;
}

int gemmTest(int output)
{
	int res = 1;
	myreal b[2][2] = {0, 0, 0, 0};
	myreal a[3][2] = {1, 2, 3, 4, 5, 6};
	int m = 2, n = 2, k = 3;
	cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, m, n, k,
		1, a, n, a, n, 0, b, n);

	res = res && (b[0][0] == 35);
	res = res && (b[0][1] == 44);
	res = res && (b[1][0] == 44);
	res = res && (b[1][1] == 56);

	if (output) {
		if (res) {
			printf("gemmTest: PASSED\n");
		} else {
			printf("gemmTest: FAILED\n");
		}
	}
	return res;
}


int gesvTest(int output)
{
	int res = 1;
	myreal a[3][3] = {1, 4, 7, 2, 5, 8, 3, 6, 0};
	myreal b[3] = {10, 11, 12};
	int n = 3;
	int nrhs = 1;
	int ipiv[3];
	int info;
	myreal prec = 1e-15;

	dgesv(&n, &nrhs, a, &n, ipiv, b, &n, &info);

	res = res && (abs(b[0] + 28.0/3) < prec);
	res = res && (abs(b[1] - 29.0/3) < prec);
	res = res && (abs(b[2]) < prec);

	if (output) {
		if (res) {
			printf("gesvTest: PASSED\n");
		} else {
			printf("gesvTest: FAILED\n");
		}
	}
	return res;
}

int cbrtTest(int output)
{
	int res = abs(cbrt(729.0, 1e-20) - 9) < 1e20;
	if (output) {
		if (res) {
			printf("cbrtTest: PASSED\n");
		} else {
			printf("cbrtTest: FAILED\n");
		}
	}
	return res;
}

int equTest(int output)
{
	int res = 1;
	res = res && solve2Test(output);
	res = res && solve3Test(output);
	res = res && solve4Test(output);
	if (output) {
		if (res) {
			printf("equTest: ALL TESTS PASSED\n");
		} else {
			printf("equTest: SOME TESTS FAILED\n");
		}
	}
	return res;
}

int solve2Test(int output)
{
	int res = 1;
	myreal prec = 1e-12;
	mycomplex tmp; 
	mycomplex a = makeComplex(1, 0);
	mycomplex b = makeComplex(5, 0);
	mycomplex c = makeComplex(6, 0);
	mycomplex r[2];
	int i;

	solve2(a, b, c, r);

	for (i = 0; i < 2; i++) {
		tmp = makeComplex(0, 0);
		tmp = zaddz(tmp, zmulz(a, zmulz(r[i], r[i]))); 
		tmp = zaddz(tmp, zmulz(b, r[i])); 
		tmp = zaddz(tmp, c);
		res = res && (zabs(tmp) < prec);
	}
	
	if (output) {
		if (res) {
			printf("solve2Test: PASSED\n");
		} else {
			printf("solve2Test: FAILED\n");
		}
	}
	return res;
}

int solve3Test(int output)
{
	int res = 1;
	myreal prec = 1e-12;
	mycomplex tmp; 
	myreal a = 1;
	myreal b = 5;
	myreal c = 6;
	myreal d = -3;
	mycomplex r[3];
	int i;

	solve3(a, b, c, d, r);

	for (i = 0; i < 3; i++) {
		tmp = makeComplex(0, 0);
		tmp = zaddz(tmp, zmulr(zmulz(r[i], zmulz(r[i], r[i])), a)); 
		tmp = zaddz(tmp, zmulr(zmulz(r[i], r[i]), b)); 
		tmp = zaddz(tmp, zmulr(r[i], c)); 
		tmp = zaddr(tmp, d);
		res = res && (zabs(tmp) < prec);
	}
	
	if (output) {
		if (res) {
			printf("solve3Test: PASSED\n");
		} else {
			printf("solve3Test: FAILED\n");
		}
	}
	return res;
}

int solve4Test(int output)
{
	int res = 1;
	myreal prec = 1e-12;
	mycomplex tmp; 
	myreal a = 1;
	myreal b = 5;
	myreal c = 6;
	myreal d = -3;
	myreal e = -21;
	mycomplex r[4];
	int i;

	solve4(a, b, c, d, e, r);

	for (i = 0; i < 4; i++) {
		tmp = makeComplex(0, 0);
		tmp = zaddz(tmp, zmulr(zmulz(zmulz(r[i], r[i]), zmulz(r[i], r[i])), a)); 
		tmp = zaddz(tmp, zmulr(zmulz(r[i], zmulz(r[i], r[i])), b)); 
		tmp = zaddz(tmp, zmulr(zmulz(r[i], r[i]), c)); 
		tmp = zaddz(tmp, zmulr(r[i], d)); 
		tmp = zaddr(tmp, e);
		res = res && (zabs(tmp) < prec);
	}
	
	if (output) {
		if (res) {
			printf("solve4Test: PASSED\n");
		} else {
			printf("solve4Test: FAILED\n");
		}
	}
	return res;
}

int distToEllipseTest(int output)
{
	int res = 1;
	myreal prec = 1e-10;
	myreal X = 3;
	myreal Y = 1;
	myreal a = 4;
	myreal b = 3;
	myreal d[2];

	distToEllipse(a, b, X, Y, d);
	res = res && (abs(d[0] - 3.5510185489) < prec);
	res = res && (abs(d[1] - 1.3809508813) < prec);

	if (output) {
		if (res) {
			printf("distToEllipseTest: PASSED\n");
		} else {
			printf("distToEllipseTest: FAILED\n");
		}
	}
	return res;
}

int ellipseFittingTest(int output)
{	
	int res = 1;
	myreal Xc = -7;
	myreal Yc = 21; 
	myreal a = 130;
	myreal b = 90;
	myreal al = M_PI / 8;
	myreal step = M_PI / 1239;
	int size = 2*M_PI / step + 1;
	myreal t;
	int i;
	myreal **X = (myreal **)malloc(size * sizeof(myreal *));
	myreal ellipseParams[5];
	myreal prec = 0.001;
	for (i = 0, t = 0; i < size; i++, t+= step) {
		myreal x = a*cos(t);
		myreal y = b*sin(t);
		X[i] = (myreal *)malloc(2 * sizeof(myreal));
		X[i][0] = cos(al)*x - sin(al)*y + Xc;
		X[i][1] = sin(al)*x + cos(al)*y + Yc;
	}
	ellipseFitting(X, size, ellipseParams);
	
	res = res && (abs(ellipseParams[0] - Xc) < prec);
	res = res && (abs(ellipseParams[1] - Yc) < prec);
	res = res && (abs(ellipseParams[2] - a) < prec);
	res = res && (abs(ellipseParams[3] - b) < prec);
	res = res && (abs(ellipseParams[4] - al) < prec);
	for (i = 0; i < size; i++) {
		free(X[i]);
	}
	free(X);
	if (output) {
		if (res) {
			printf("ellipseFittingTest: PASSED\n");
		} else {
			printf("ellipseFittingTest: FAILED\n");
		}
	}
	return res;
}

int ellipseFittingCircleTest(int output)
{
	int res = 1;
	myreal Xc = 0;
	myreal Yc = 0;
	myreal R = 10;
	myreal step = M_PI / 1239;
	int size = 2 * M_PI / step + 1;
	myreal t;
	int i;
	myreal **X = (myreal **)malloc(size * sizeof(myreal *));
	myreal ellipseParams[5];
	myreal prec = 0.001;
	for (i = 0, t = 0; i < size; i++, t += step) {
		myreal x = R*cos(t);
		myreal y = R*sin(t);
		X[i] = (myreal *)malloc(2 * sizeof(myreal));
		X[i][0] = x + Xc;
		X[i][1] = y + Yc;
	}
	ellipseFitting(X, size, ellipseParams);

	res = res && (abs(ellipseParams[0] - Xc) < prec);
	res = res && (abs(ellipseParams[1] - Yc) < prec);
	res = res && (abs(ellipseParams[2] - R) < prec);
	res = res && (abs(ellipseParams[3] - R) < prec);
	res = res && (abs(ellipseParams[4] - 0) < prec);
	for (i = 0; i < size; i++) {
		free(X[i]);
	}
	free(X);
	if (output) {
		if (res) {
			printf("ellipseFittingCircleTest: PASSED\n");
		}
		else {
			printf("ellipseFittingCircleTest: FAILED\n");
		}
	}
	return res;
}

void ellipseFittingNoiseTest(int output)
{
	int res = 1;
	myreal Xc = -7;
	myreal Yc = 21;
	myreal a = 130;
	myreal b = 90;
	myreal al = -M_PI / 6;
	myreal step = M_PI / 120;
	int size = 2 * M_PI / step + 1;
	myreal t;
	int i;
	myreal **X = (myreal **)malloc(size * sizeof(myreal *));
	myreal ellipseParams[5];
	myreal prec = 0.001;
	myreal noise_level = 0.1;
	srand(time(0));
	for (i = 0, t = 0; i < size; i++, t += step) {
		myreal x = a*cos(t);
		myreal y = b*sin(t);
		X[i] = (myreal *)malloc(2 * sizeof(myreal));
		myreal r0 = noise_level * (myreal)(RAND_MAX/2 - rand()) / RAND_MAX;
		myreal r1 = noise_level * (myreal)(RAND_MAX/2 - rand()) / RAND_MAX;
		X[i][0] = cos(al)*x - sin(al)*y + Xc;
		X[i][0] += X[i][0] * r0;
		X[i][1] = sin(al)*x + cos(al)*y + Yc;
		X[i][1] += X[i][1] * r1;
	}
	ellipseFitting(X, size, ellipseParams);

	res = res && (abs(ellipseParams[0] - Xc) < prec);
	res = res && (abs(ellipseParams[1] - Yc) < prec);
	res = res && (abs(ellipseParams[2] - a) < prec);
	res = res && (abs(ellipseParams[3] - b) < prec);
	res = res && (abs(ellipseParams[4] - al) < prec);
	for (i = 0; i < size; i++) {
		free(X[i]);
	}
	free(X);
	if (output) {
		if (res) {
			printf("ellipseFittingNoiseTest: PASSED\n");
		}
		else {
			printf("ellipseFittingNoiseTest: FAILED\n");
		}
	}
	return res;
}