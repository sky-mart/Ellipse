#ifndef MYCOMPLEX_H
#define MYCOMPLEX_H

#include <math.h>

#define M_PI       3.14159265358979323846
#define M_PI_2     1.57079632679489661923
#define M_PI_4     0.785398163397448309616

typedef double myreal;
typedef struct {
	myreal r;
	myreal i;
} mycomplex;

mycomplex makeComplex(myreal r, myreal i);

mycomplex conj(mycomplex z);
myreal zabs(mycomplex z);
myreal zabs2(mycomplex z);
myreal sign(myreal n);

mycomplex toPolar(mycomplex cart);
mycomplex toCart(mycomplex polar);

mycomplex complexNthroot(mycomplex z, int n);

mycomplex zaddz(mycomplex z1, mycomplex z2);
mycomplex zaddr(mycomplex z, myreal r);
mycomplex zsubz(mycomplex z1, mycomplex z2);
mycomplex zsubr(mycomplex z, myreal r);
mycomplex zmulz(mycomplex z1, mycomplex z2);
mycomplex zmulr(mycomplex z, myreal r);
mycomplex zdivz(mycomplex z1, mycomplex z2);
mycomplex zdivr(mycomplex z, myreal r);

void printComplex(mycomplex z);

#endif