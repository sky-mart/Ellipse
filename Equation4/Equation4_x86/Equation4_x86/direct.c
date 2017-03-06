//
//  direct.c
//  Equation4_x86
//
//  Created by Влад Агиевич on 06.03.17.
//  Copyright © 2017 Mocsmart. All rights reserved.
//

#include "direct.h"

#define SQRT3 1.732050807568877f

/*
 * solves equation ax^2 + bx + c = 0
 * puts roots in the array
 * returns 1 on success
 */
int direct_solve2(dreal a, dreal b, dreal c, dcomplex roots[2])
{
    dreal D;
    if (a == 0) {
        if (b == 0) {
            return 0; // any x is a solution
        }
        roots[0] = roots[1] = -c / b;
        return 1;
    }
    
    D = b*b - 4*a*c;
    if (D >= 0) {
        roots[0] = (-b + sqrtd(D)) / (2*a);
        roots[1] = (-b - sqrtd(D)) / (2*a);
    } else {
        roots[0] = (-b + sqrtd(-D) * I) / (2*a);
        roots[1] = (-b - sqrtd(-D) * I) / (2*a);
    }
    return 1;
}

/*
 * solves equation ax^3 + bx^2 + cx + d = 0
 * puts roots in the array
 * returns 1 on success
 */
int direct_solve3(dreal a, dreal b, dreal c, dreal d, dcomplex roots[3])
{
    dreal p, q, Q;
    dcomplex alpha, beta;
    
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
        alpha   = cbrtd(-q/2 + sqrtd(Q));
        beta    = cbrtd(-q/2 - sqrtd(Q));
    } else {
        alpha   = cpowd(-q/2 + sqrtd(-Q) * I, 1.0/3);
        beta    = cpowd(-q/2 - sqrtd(-Q) * I, 1.0/3);
    }
    
    roots[0] = alpha + beta - b/(3*a);
    roots[1] = -(alpha+beta)/2 - b/(3*a) + (alpha-beta)*SQRT3/2 * I;
    roots[2] = -(alpha+beta)/2 - b/(3*a) - (alpha-beta)*SQRT3/2 * I;
    
    return 1;
}

/*
 * solves equation ax^4 + bx^3 + cx^2 + dx + e = 0
 * puts roots in the array
 * returns 1 on success
 */
int direct_solve4(dreal a, dreal b, dreal c, dreal d, dreal e, dcomplex roots[4])
{
    dreal p, q, r;
    dreal A, B, C, D, s;
    dcomplex cube[3];
    dreal a1, b1, c1, a2, b2, c2;
    int i;
    
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
    if (!direct_solve3(A, B, C, D, cube)) {
        return 0;
    }
    
    s = 0;
    if (creald(cube[0]) > p/2) {
        s = creald(cube[0]);
    } else if (creald(cube[1]) > p/2) {
        s = creald(cube[1]);
    } else if (creald(cube[2]) > p/2) {
        s = creald(cube[2]);
    }

    a1 = 1; b1 = -sqrtd(2*s-p); c1 = q/(2*sqrtd(2*s-p)) + s;
    a2 = 1; b2 = sqrtd(2*s-p); c2 = -q/(2*sqrtd(2*s-p)) + s;
    
    if (!direct_solve2(a1, b1, c1, &roots[0]) || !direct_solve2(a2, b2, c2, &roots[2])) {
        return 0;
    }
    for (i = 0; i < 4; i++) {
        roots[i] -= a/4;
    }

    return 1;
}

