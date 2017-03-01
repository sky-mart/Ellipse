//
//  equations.c
//  Eq4SinglePrec
//
//  Created by ???? ??????? on 28.02.17.
//  Copyright © 2017 Mocsmart. All rights reserved.
//

#include "equations.h"
#include "extra_math.h"

/*
 * solves equation ax^2 + bx + c = 0
 * puts roots in the array
 * returns 1 on success
 */
int solve2(real a, real b, real c, complexr roots[2])
{
    real D;
    if (a == 0) {
        if (b == 0) {
            return 0; // any x is a solution
        }
        roots[0] = roots[1] = -c / b;
        return 1;
    }
    
    D = b*b - 4*a*c;
    if (D >= 0) {
        roots[0] = (-b + sqrtr(D)) / (2*a);
        roots[1] = (-b - sqrtr(D)) / (2*a);
    } else {
        roots[0] = (-b + sqrtr(-D) * I) / (2*a);
        roots[1] = (-b - sqrtr(-D) * I) / (2*a);
    }
    return 1;
}

/*
 * solves equation ax^3 + bx^2 + cx + d = 0
 * puts roots in the array
 * returns 1 on success
 */
int solve3(real a, real b, real c, real d, complexr roots[3])
{
    double p, q, Q;
    double complex alpha, beta;
    //    real rp[5], ip[5]; // debug
#ifdef SINGLE_PREC_DEBUG
    printf("solve3:\n");
    printf("a=%e\n", a);
    printf("b=%e\n", b);
    printf("c=%e\n", c);
    printf("d=%e\n", d);
    printf("\n");
#endif
    // reduce to the form y^3 + py + q = 0
    // substitution: x = y - b/3a
    p = c/a - b*b/(3*a*a);
    q = 2*b*b*b/(27*a*a*a) - b*c/(3*a*a) + d/a;
    
    Q = p*p*p/27 + q*q/4;
#ifdef SINGLE_PREC_DEBUG
    printf("p=%e\n", p);
    printf("q=%e\n", q);
    printf("Q=%e\n", Q);
    printf("\n");
#endif
    // Q > 0 - one real root and two complex conjugated roots
    // Q = 0 - one single real root and one double real root, or,
    //         if p = q = 0, then one triple real root
    // Q < 0 - three real roots
    
    if (Q >= 0) {
        //alpha   = mypowr(-q/2 + sqrtr(Q), 1.0/3);
        //beta    = mypowr(-q/2 - sqrtr(Q), 1.0/3);
        alpha   = cbrt(-q/2 + sqrt(Q));
        beta    = cbrt(-q/2 - sqrt(Q));
    } else {
        alpha   = cpow(-q/2 + sqrt(-Q) * I, 1.0/3);
        beta    = cpow(-q/2 - sqrt(-Q) * I, 1.0/3);
    }
#ifdef SINGLE_PREC_DEBUG
    printf("alpha="); cprint(alpha); printf("\n");
    printf("beta="); cprint(beta); printf("\n");
    printf("\n");
#endif
    //    cparts(alpha, &rp[0], &ip[0]);
    //    cparts(beta, &rp[1], &ip[1]);
    
    roots[0] = alpha + beta - b/(3*a);
    roots[1] = -(alpha+beta)/2 - b/(3*a) + (alpha-beta)*SQRT3/2 * I;
    roots[2] = -(alpha+beta)/2 - b/(3*a) - (alpha-beta)*SQRT3/2 * I;
#ifdef SINGLE_PREC_DEBUG
    printf("roots[0]="); cprint(roots[0]); printf("\n");
    printf("roots[1]="); cprint(roots[1]); printf("\n");
    printf("roots[2]="); cprint(roots[2]); printf("\n");
    printf("\n\n");
#endif
    
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
int solve4(real a, real b, real c, real d, real e, complexr roots[4])
{
    real p, q, r;
    real A, B, C, D, s;
    complexr cube[3];
    real a1, b1, c1, a2, b2, c2;
    int i;
    //real rp[3], ip[3]; // debug
#ifdef SINGLE_PREC_DEBUG
    printf("solve4:\n");
    printf("a=%e\n", a);
    printf("b=%e\n", b);
    printf("c=%e\n", c);
    printf("d=%e\n", d);
    printf("e=%e\n", e);
    printf("\n");
#endif
    
    b /= a; c /= a; d /= a; e /= a;
    a = b; b = c; c = d; d = e;
    
    // reduce to the form y^4 + p*y^2 + q*y + r = 0
    p = b - 3*a*a/8;
    q = a*a*a/8 - a*b/2 + c;
    r = - 3*a*a*a*a/256 + a*a*b/16 - c*a/4 + d;
#ifdef SINGLE_PREC_DEBUG
    printf("p=%e\n", p);
    printf("q=%e\n", q);
    printf("r=%e\n", r);
    printf("\n");
#endif
    // obtain cubic resolvent A*s^3 + B*s^2 + C*s + D = 0
    A = 2;
    B = -p;
    C = -2*r;
    D = r*p - q*q/4;
    if (!solve3(A, B, C, D, cube)) {
        return 0;
    }
    
    //    for (i = 0; i < 3; i++) {
    //        rp[i] = crealr(cube[i]);
    //        ip[i] = cimagr(cube[i]);
    //    }
    
    s = 0;
    if (crealr(cube[0]) > p/2) {
        s = crealr(cube[0]);
    } else if (crealr(cube[1]) > p/2) {
        s = crealr(cube[1]);
    } else if (crealr(cube[2]) > p/2) {
        s = crealr(cube[2]);
    }
#ifdef SINGLE_PREC_DEBUG
    printf("s=%e\n", s);
    printf("\n");
#endif
    a1 = 1; b1 = -sqrtr(2*s-p); c1 = q/(2*sqrtr(2*s-p)) + s;
    a2 = 1; b2 = sqrtr(2*s-p); c2 = -q/(2*sqrtr(2*s-p)) + s;
#ifdef SINGLE_PREC_DEBUG
    printf("a1=%e\n", a1);
    printf("b1=%e\n", b1);
    printf("c1=%e\n", c1);
    printf("a2=%e\n", a2);
    printf("b2=%e\n", b2);
    printf("c2=%e\n", c2);
    printf("\n");
#endif
    
    if (!solve2(a1, b1, c1, &roots[0]) || !solve2(a2, b2, c2, &roots[2])) {
        return 0;
    }
    for (i = 0; i < 4; i++) {
        roots[i] -= a/4;
    }
#ifdef SINGLE_PREC_DEBUG
    printf("roots[0]="); cprint(roots[0]); printf("\n");
    printf("roots[1]="); cprint(roots[1]); printf("\n");
    printf("roots[2]="); cprint(roots[2]); printf("\n");
    printf("roots[3]="); cprint(roots[3]); printf("\n");
    printf("\n\n");
#endif
    return 1;
}
