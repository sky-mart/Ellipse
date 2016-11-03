#include "mycomplex.h"
#include <stdio.h>
#include <math.h>

mycomplex makeComplex(myreal r, myreal i)
{
	mycomplex z;
	z.r = r;
	z.i = i;
	return z;
}

mycomplex conj(mycomplex z)
{
	mycomplex c;
	c.r = z.r;
	c.i = -z.i;
	return c;
}

myreal zabs(mycomplex z)
{
	mycomplex a = zmulz(z, conj(z));
	return sqrt(a.r);
}

myreal zabs2(mycomplex z)
{
	mycomplex a = zmulz(z, conj(z));
	return a.r;
}

mycomplex complexNthroot(mycomplex z, int n)
{
	mycomplex polar = toPolar(z);
	polar.r = exp(log(polar.r) / n);
	polar.i = polar.i / n;
	return toCart(polar);
}

myreal sign(myreal n)
{
	if (n > 0) {
		return 1;
	} if (n < 0) {
		return -1;
	} return 0;
}
mycomplex toPolar(mycomplex cart)
{
	mycomplex polar;
	double signX = sign(cart.r);
	double signY = sign(cart.i);
	polar.r = sqrt(cart.r*cart.r + cart.i*cart.i);
	if (signY) {
		polar.i = atan((cart.i) / (cart.r)) + M_PI_2 * signY * (1 - signX);
	} else {
		polar.i = signX >= 0 ? 0 : M_PI; 
	}
	return polar;
}

mycomplex toCart(mycomplex polar)
{
	mycomplex cart;
	cart.r = polar.r*cos(polar.i);
	cart.i = polar.r*sin(polar.i);
	return cart;
}

mycomplex zaddz(mycomplex z1, mycomplex z2)
{
	return makeComplex(z1.r + z2.r, z1.i + z2.i);
}

mycomplex zaddr(mycomplex z, myreal r)
{
	return makeComplex(z.r + r, z.i);
}

mycomplex zsubz(mycomplex z1, mycomplex z2)
{
	return makeComplex(z1.r - z2.r, z1.i - z2.i);
}

mycomplex zsubr(mycomplex z, myreal r)
{
	return makeComplex(z.r - r, z.i);
}

mycomplex zmulz(mycomplex z1, mycomplex z2)
{
	return makeComplex(z1.r*z2.r - z1.i*z2.i, z1.r*z2.i + z1.i*z2.r);
}

mycomplex zmulr(mycomplex z, myreal r)
{
	return makeComplex(z.r*r, z.i*r);
}

mycomplex zdivz(mycomplex z1, mycomplex z2)
{
	return zdivr(zmulz(z1, conj(z2)), zabs2(z2));
}

mycomplex zdivr(mycomplex z, myreal r)
{
	return makeComplex(z.r/r, z.i/r);
}

void printComplex(mycomplex z)
{
	printf("%lf", z.r);
	if (z.i >= 0) {
		printf("+%lfi", z.i);
	} else {
		printf("%lfi", z.i);
	}
}