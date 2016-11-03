#ifndef ELLIPSE_H
#define ELLIPSE_H

#include "mycomplex.h"

void ellipseFitting(myreal **X, int n, myreal ellipseParams[5]);

void jacobianAndPenalty(myreal Xi[2], myreal ellipseParams[5], myreal Ji[2][5], myreal ei[2]);

void linsolve(myreal *J, myreal *e, myreal *x, int nRows, int nParams);

void massCenter(myreal **X, int n, myreal params[3]);

void conversion(myreal X[2], myreal T[2][2], myreal Xc[2], myreal x[2]);
void invConversion(myreal x[2], myreal T[2][2], myreal Xc[2], myreal X[2]);

void B(myreal a, myreal b, myreal C, myreal S, myreal xs[2], myreal x[2], myreal M[2][5]);

void simpleInv(myreal A[2][2], myreal invA[2][2]);

void distToEllipse(myreal a, myreal b, myreal X0, myreal Y0, myreal d[2]);

void solve4(myreal a0, myreal a1, myreal a2, myreal a3, myreal a4, mycomplex roots[4]);

void solve3(myreal a0, myreal a1, myreal a2, myreal a3, mycomplex roots[3]);

void solve2(mycomplex a0, mycomplex a1, mycomplex a2, mycomplex roots[2]);

//myreal cbrt(myreal a, myreal e);

#endif