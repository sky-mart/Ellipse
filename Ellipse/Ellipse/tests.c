#include <stdio.h>

#include "tests.h"
#include "mycomplex.h"
#include "ellipse.h"
#include "mkl.h"
#include <stdlib.h>
#include <time.h>

myreal angleTo0_2pi(myreal x, myreal prec){
	while (x < 0 && fabs(x) > prec) {
		x += 2 * M_PI;
	}
	while (x >= 2 * M_PI && fabs(x - 2 * M_PI) > prec) {
		x -= 2 * M_PI;
	}
	return x;
}

int compareAngles(myreal alpha, myreal beta, myreal prec)
{
	alpha = angleTo0_2pi(alpha, prec);
	beta = angleTo0_2pi(beta, prec);
	return fabs(alpha - beta) < prec;
}

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

void ellipseFittingPoints(myreal **X, int points_num, myreal Xc, myreal Yc, myreal a, myreal b, myreal alpha, myreal noise_level)
{	
	myreal step = 2 * M_PI / (points_num - 1);
	myreal t, x, y;
	myreal rand_x, rand_y;
	int i;
	srand(time(0));
	for (i = 0, t = 0; i < points_num; i++, t += step) {
		x = a*cos(t);
		y = b*sin(t);
		rand_x = noise_level * (myreal)(RAND_MAX / 2 - rand()) / RAND_MAX;
		rand_y = noise_level * (myreal)(RAND_MAX / 2 - rand()) / RAND_MAX;
		X[i][0] = (cos(alpha)*x - sin(alpha)*y + Xc) * (1 + rand_x);
		X[i][1] = (sin(alpha)*x + cos(alpha)*y + Yc) * (1 + rand_y);
	}
}

int ellipseFittingCompare(myreal ellipseParams[5], myreal Xc, myreal Yc, myreal a, myreal b, myreal alpha, myreal prec)
{
	int res = 1;
	res = res && (abs(ellipseParams[0] - Xc) < prec);
	res = res && (abs(ellipseParams[1] - Yc) < prec);

	if ((fabs(ellipseParams[2] - a) < prec) && (fabs(ellipseParams[3] - b) < prec)) {
		if (fabs(a - b) > prec)
			res = res && (compareAngles(alpha, ellipseParams[4], prec) |
			compareAngles(alpha + M_PI, ellipseParams[4], prec));
	}
	else if ((fabs(ellipseParams[2] - b) < prec) && (fabs(ellipseParams[3] - a) < prec)) {
		if (fabs(a - b) > prec)
			res = res && (compareAngles(alpha + M_PI_2, ellipseParams[4], prec) |
			compareAngles(alpha + 3 * M_PI_2, ellipseParams[4], prec));
	}
	else {
		res = 0;
	}
	return res;
}

int ellipseFittingTest(myreal **X, int points_num, myreal Xc, myreal Yc, myreal a, myreal b, myreal alpha, myreal noise_level, myreal prec)
{
	int res = 1;
	myreal ellipseParams[5];
	ellipseFittingPoints(X, points_num, Xc, Yc, a, b, alpha, noise_level);
	ellipseFitting(X, points_num, ellipseParams);

	return ellipseFittingCompare(ellipseParams, Xc, Yc, a, b, alpha, prec);
}

void ellipseFittingOutput(int res, int points_num, myreal Xc, myreal Yc, myreal a, myreal b, myreal alpha, myreal noise_level, myreal prec)
{
	if (res) {
		printf("Ellipse Fitting(N = %d, noise= %.2lf, prec= %.2e, Xc = %.2lf, Yc = %.2lf, a = %.2lf, b = %.2lf, alpha = %.2lf) : PASSED\n",
			points_num, noise_level, prec, Xc, Yc, a, b, alpha);
	}
	else {
		printf("Ellipse Fitting(N = %d, noise= %.2lf, prec= %.2e, Xc = %.2lf, Yc = %.2lf, a = %.2lf, b = %.2lf, alpha = %.2lf) : FAILED\n",
			points_num, noise_level, prec, Xc, Yc, a, b, alpha);
	}
}

int ellipseFittingTestWrapper(int points_num, myreal Xc, myreal Yc, myreal a, myreal b, myreal alpha, myreal noise_level, myreal prec, int output)
{
	int i, res;
	myreal **X = (myreal **)malloc(points_num * sizeof(myreal *));
	for (i = 0; i < points_num; i++) {
		X[i] = (myreal *)malloc(2 * sizeof(myreal));
	}
	res = ellipseFittingTest(X, points_num, Xc, Yc, a, b, alpha, noise_level, prec);
	if (output)
		ellipseFittingOutput(res, points_num, Xc, Yc, a, b, alpha, noise_level, prec);

	for (i = 0; i < points_num; i++) {
		free(X[i]);
	}
	free(X);
}

int ellipseFittingLoopTest(int output)
{
	int res;
	myreal Xc, Yc, a, b, alpha;
	int points_num;
	int i;
	myreal ellipseParams[5];
	myreal prec = 0.001;
	myreal **X;

	a = 1.0;
	for (points_num = 30; points_num <= 1000; points_num += 1) {
		X = (myreal **)malloc(points_num * sizeof(myreal *));
		for (i = 0; i < points_num; i++) {
			X[i] = (myreal *)malloc(2 * sizeof(myreal));
		}
		
		for (Xc = -1.0; Xc <= 1.0; Xc += 0.1) {
			for (Yc = -1.0; Yc <= 1.0; Yc += 0.1) {
				for (b = 1.0; b <= 10.0; b += 0.5) {
					for (alpha = 0.0; alpha <= M_PI_2; alpha += 0.2) {
						myreal noise_level = 0;
						ellipseFittingPoints(X, points_num, Xc, Yc, a, b, alpha, noise_level);
						res = ellipseFittingTest(X, points_num, Xc, Yc, a, b, alpha, noise_level, prec);

						if (output) ellipseFittingOutput(res, points_num, Xc, Yc, a, b, alpha, noise_level, prec);

						if (!res) goto quit;
					}
				}
			}
		}
		for (i = 0; i < points_num; i++) {
			free(X[i]);
		}
		free(X);
	}
quit:
	printf("Test finished\n");
}