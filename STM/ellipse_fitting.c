#include "ellipse_fitting.h"

#include "math.h"
#include "float.h"

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
