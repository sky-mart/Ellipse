//
//  fourier_sturm.h
//  Equation4_x86
//
//  Created by Влад Агиевич on 07.03.17.
//  Copyright © 2017 Mocsmart. All rights reserved.
//

#ifndef fourier_sturm_h
#define fourier_sturm_h

#include "settings.h"

typedef struct {
    real from;
    real to;
} segment;

typedef struct {
    real *a;
    uint n;
} poly;

void segment_print(segment *s);
void segment_println(segment *s);


int     poly_init(poly *p, uint n);
void    poly_free(poly *p);
real    poly_value(poly *p, real x);
void    poly_root_segments(poly *p, segment *neg, segment *pos);
int     poly_div(poly *p, poly *d, poly *q, poly *r);
void    poly_diff(poly *p, poly *d);
void    poly_print(poly *p);
void    poly_println(poly *p);
int     poly_sturmseq(poly *p, poly *seq);
uint poly_sturmseq_signchangenum(poly *seq, real x);
int poly_sturm_theorem(poly *p, segment *s, uint *root_count);
void poly_sturm_using_seq(poly *seq, segment *s, uint *root_count);
int poly_localize(poly *p, segment *root_segments, uint *root_count);
int poly_solve(poly *p, real *roots, uint *root_count);

int sturm_solve4(real a, real b, real c, real d, real e, real *roots, uint *root_count);
int real_root_count4(real a, real b, real c, real d, real e, uint *root_count);

#endif /* fourier_sturm_h */
