#ifndef TESTS_H
#define TESTS_H

int complexTest(int output);
int addTest(int output);

int blasTest(int output);
int copyTest(int output);
int axpyTest(int output);
int gemvTest(int output);
int gemmTest(int output);

int gesvTest(int output);

int cbrtTest(int output);

int equTest(int output);
int solve2Test(int output);
int solve3Test(int output);
int solve4Test(int output);

int distToEllipseTest(int output);
int ellipseFittingTest(int output);
int ellipseFittingCircleTest(int output);

#endif