//
//  fourier_sturm.c
//  Equation4_x86
//
//  Created by Влад Агиевич on 07.03.17.
//  Copyright © 2017 Mocsmart. All rights reserved.
//

#include "sturm.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

real segment_len(segment *s)
{
    return s->to - s->from;
}

void segment_print(segment *s)
{
    printf("(%f, %f)", s->from, s->to);
}

void segment_println(segment *s)
{
    segment_print(s);
    printf("\n");
}

void ring_buffer_incr(uint size, uint *pIndex)
{
    (*pIndex)++;
    if (*pIndex == size)
        *pIndex = 0;
}

uint ring_buffer_len(uint size, uint begin_index, uint end_index)
{
    if (end_index > begin_index)
        return end_index - begin_index;
    return end_index + size - begin_index;
}

int ring_buffer_normalize(segment **pBuf, uint size, uint begin_index, uint end_index)
{
    segment *tmp_buf;
    uint i, len;
    if ( !(tmp_buf = (segment *)malloc(size * sizeof(segment))) )
        return 0;
    
    if (end_index > begin_index) {
        len = end_index - begin_index;
        memcpy(tmp_buf, &(*pBuf)[begin_index], len * sizeof(segment));
    } else {
        memcpy(tmp_buf, &(*pBuf)[begin_index], (size - begin_index) * sizeof(segment));
        memcpy(&tmp_buf[size-begin_index], &(*pBuf)[0], end_index * sizeof(segment));
        len = end_index + size - begin_index;
    }
    
    for (i = 0; i < len; i++)
        (*pBuf)[i] = tmp_buf[i];
    
    free(tmp_buf);
    return 1;
}

void poly_print(poly *p)
{
    uint i, pow;
    for (i = 0; i <= p->n; i++) {
        pow = p->n - i;
        
        if (p->a[i] != 0) {
            if (i != 0) {
                if (p->a[i] > 0)
                    printf("+");
            }
            printf("%f", p->a[i]);
            if (i != p->n) {
                printf("*x");
                if (i != p->n -1)
                    printf("^%i", pow);
            }
        }
    }
}

void poly_println(poly *p)
{
    poly_print(p);
    printf("\n");
}

void poly_sturmseq_println(poly *seq)
{
    uint i;
    for (i = 0; i <= seq[0].n; i++)
        poly_println(&seq[i]);
}

void poly_swap_cf(poly *u, poly *v)
{
    real *tmp;
    tmp = u->a;
    u->a = v->a;
    v->a = tmp;
}

int poly_init(poly *p, uint n)
{
    if ( !(p->a = (real *)malloc((n + 1) * sizeof(real))) ) {
        return 0;
    }
    p->n = n;
    memset(p->a, 0, (n + 1) * sizeof(real));
    return 1;
}

void poly_free(poly *p)
{
    free(p->a);
}

real poly_value(poly *p, real x)
{
    int i;
    real val = 0, xpow = 1;
    for (i = p->n; i >= 0; i--) {
        val += p->a[i] * xpow;
        xpow *= x;
    }
    return val;
}

// d = p'
void poly_diff(poly *p, poly *d)
{
    uint i;
    d->n = p->n - 1;
    for (i = 0; i <= d->n; i++) {
        d->a[i] = p->a[i] * (p->n - i);
    }
}

// just copies polynomial's coefficients and power
void poly_copy(poly *src, poly *dst)
{
    dst->n = src->n;
    memcpy(dst->a, src->a, (src->n + 1) * sizeof(real));
}

// as for now polynomial powers are supposed equal
void poly_sub(poly *p, poly *subtrahend, poly *subres)
{
    uint i, j;
    for (i = 0; i <= p->n; i++) {
        subres->a[i] = p->a[i] - subtrahend->a[i];
    }
    for (i = 0; i <= p->n; i++) {
        if (subres->a[i] != 0) {
            if (i != 0) {
                for (j = i; j <= p->n; j++)
                    subres->a[j-i] = subres->a[j];
                subres->n -= i;
            }
            break;
        }
    }
}

void poly_mulc(poly *p, real c, poly *res)
{
    uint i;
    res->n = p->n;
    for (i = 0; i <= p->n; i++) {
        res->a[i] = p->a[i] * c;
    }
}

// p = q*d + r
int poly_div(poly *p, poly *d, poly *q, poly *r)
{
    uint i;
    poly subtrahend, subres;
    real c;
    
    if (!poly_init(&subtrahend, p->n) || !poly_init(&subres, p->n))
        return 0;
    
    poly_copy(p, r);
    q->n = p->n - d->n;
    
    while (r->n >= d->n) {
        i = p->n - r->n;
        subtrahend.n = r->n;
        c = r->a[0] / d->a[0];
        poly_mulc(d, c, &subtrahend);
        poly_sub(r, &subtrahend, &subres);
        
        q->a[i] = c;
        poly_swap_cf(r, &subres);
        r->n = subres.n;
    }
    
    poly_free(&subtrahend);
    poly_free(&subres);
    return 1;
}

int poly_sturmseq(poly *p, poly *seq)
{
    poly q, r;
    uint i;
    
    if (!poly_init(&q, 1) || !poly_init(&r, p->n))
        return 0;

    poly_copy(p, &seq[0]);
    poly_diff(p, &seq[1]);
    
    for (i = 1; i < p->n; i++) {
        if (!poly_div(&seq[i-1], &seq[i], &q, &r))
            return 0;
        poly_mulc(&r, -1, &seq[i+1]);
    }
    
    poly_free(&q);
    poly_free(&r);
    return 1;
}

int poly_sturmseq_init(poly *p, poly **pSeq)
{
    uint i;
    *pSeq = (poly *)malloc((p->n + 1) * sizeof(poly));
    for (i = 0; i <= p->n; i++) {
        if (!poly_init(&(*pSeq)[i], p->n - i))
            return 0;
    }
    if (!poly_sturmseq(p, *pSeq))
        return 0;
    return 1;
}

void poly_sturmseq_free(poly *seq)
{
    uint i;
    for (i = 0; i <= seq[0].n; i++) {
        poly_free(&seq[i]);
    }
    free(seq);
}

uint poly_sturmseq_signchangenum(poly *seq, real x)
{
    uint i;
    uint counter = 0;
    real val1, val2;
    
    val1 = poly_value(&seq[0], x);
    for (i = 1; i <= seq[0].n; i++) {
        val2 = poly_value(&seq[i], x);
        if (val2 != 0) {
            counter += (val1 * val2 < 0) ? 1 : 0;
            val1 = val2;
        }
    }
    return counter;
}

int poly_sturm_theorem(poly *p, segment *s, uint *root_count)
{
    poly *seq;
    
    if (!poly_sturmseq_init(p, &seq))
        return 0;
    poly_sturmseq_println(seq);
    poly_sturm_using_seq(seq, s, root_count);
    
    poly_sturmseq_free(seq);
    return 1;
}

void poly_sturm_using_seq(poly *seq, segment *s, uint *root_count)
{
    uint from_count, to_count;
    from_count  = poly_sturmseq_signchangenum(seq, s->from);
    to_count    = poly_sturmseq_signchangenum(seq, s->to);
    *root_count  = from_count - to_count;
}

void poly_root_segments(poly *p, segment *neg, segment *pos)
{
    uint i;
    real big, tmp, abs_a0, abs_an;
    real A, B;
    
    // maximum of a[1] - a[n-1]
    big = fabsr(p->a[1]);
    for (i = 2; i < p->n; i++) {
        if ( (tmp = fabsr(p->a[i])) > big)  {
            big = tmp;
        }
    }
    abs_a0 = fabsr(p->a[0]);
    abs_an = fabsr(p->a[p->n]);
    A = abs_an > big ? abs_an : big;
    B = abs_a0 > big ? abs_a0 : big;
    
    pos->from = abs_an / (abs_an + B);
    pos->to = 1 + A / abs_a0;
    
    neg->from = - pos->to;
    neg->to = - pos->from;
}

int poly_localize_roots(poly *p, segment *root_segments, uint *root_count)
{
    poly *seq;
    segment neg, pos;
    uint seg_begin_index, seg_end_index, i;
    uint *root_count_inseg, *w_a, *w_b;
    real seg_len;
    
    if (!(root_count_inseg = (uint *)malloc(p->n * sizeof(uint))) ||
        !(w_a = (uint *)malloc(p->n * sizeof(uint))) ||
        !(w_b = (uint *)malloc(p->n * sizeof(uint))) )
            return 0;
    
    if (!poly_sturmseq_init(p, &seq))
        return 0;
    
    poly_root_segments(p, &neg, &pos);
    
    seg_begin_index = 0;
    seg_end_index = 0;
    *root_count = 0;
    
    w_a[seg_end_index] = poly_sturmseq_signchangenum(seq, neg.from);
    w_b[seg_end_index] = poly_sturmseq_signchangenum(seq, neg.to);
    root_count_inseg[seg_end_index] = w_a[seg_end_index] - w_b[seg_end_index];
    if (root_count_inseg[seg_end_index] > 0) {
        *root_count += root_count_inseg[seg_end_index];
        root_segments[seg_end_index] = neg;
        seg_end_index++;
    }
    
    w_a[seg_end_index] = poly_sturmseq_signchangenum(seq, pos.from);
    w_b[seg_end_index] = poly_sturmseq_signchangenum(seq, pos.to);
    root_count_inseg[seg_end_index] = w_a[seg_end_index] - w_b[seg_end_index];
    if (root_count_inseg[seg_end_index] > 0) {
        *root_count += root_count_inseg[seg_end_index];
        root_segments[seg_end_index] = pos;
        seg_end_index++;
    }
    
    for (i = seg_begin_index; i < seg_end_index; i++) {
        if (ring_buffer_len(p->n, seg_begin_index, seg_end_index) == *root_count)
            break;
        
        // (*pRoot_segments)[i] - cur_seg
        if (root_count_inseg[i] > 1) {
            seg_len = segment_len(&root_segments[i]);
            
            root_segments[seg_end_index] = root_segments[i];
            root_segments[seg_end_index].to  -= seg_len / 2;
            w_a[seg_end_index] = w_a[i];
            w_b[seg_end_index] = poly_sturmseq_signchangenum(seq, root_segments[seg_end_index].to);
            root_count_inseg[seg_end_index] = w_a[seg_end_index] - w_b[seg_end_index];
            
            if (root_count_inseg[seg_end_index] > 0) {
                ring_buffer_incr(p->n, &seg_end_index);
            }
            
            root_segments[seg_end_index] = root_segments[i];
            root_segments[seg_end_index].from  += seg_len / 2;
            w_a[seg_end_index] = poly_sturmseq_signchangenum(seq, root_segments[seg_end_index].from);
            w_b[seg_end_index] = w_b[i];
            root_count_inseg[seg_end_index] = w_a[seg_end_index] - w_b[seg_end_index];
            
            if (root_count_inseg[seg_end_index] > 0) {
                ring_buffer_incr(p->n, &seg_end_index);
            }
            
            ring_buffer_incr(p->n, &seg_begin_index);
        }
    }
    
    if (!ring_buffer_normalize(&root_segments, p->n, seg_begin_index, seg_end_index))
        return 0;
    
    poly_sturmseq_free(seq);
    free(root_count_inseg);
    free(w_a);
    free(w_b);
    
    return 1;
}

real poly_find_root(poly *p, segment *s, real prec)
{
    real from, mid, to;
    real val_from, val_mid, val_to;
    
    from = s->from;
    to = s->to;
    val_from = poly_value(p, from);
    val_to = poly_value(p, to);
    
    while ((to - from) > prec) {
        mid = (from + to) / 2;
        val_mid = poly_value(p, mid);
        
        if (val_from * val_mid < 0) {
            to = mid;
            val_to = val_mid;
        } else {
            from = mid;
            val_from = val_mid;
        }
    }
    return from;
}

int poly_solve(poly *p, real *roots, uint *root_count)
{
    segment *root_segments;
    uint i;
    
    root_segments = (segment *)malloc(p->n * sizeof(segment));
    poly_localize_roots(p, root_segments, root_count);
    for (i = 0; i < *root_count; i++) {
        roots[i] = poly_find_root(p, &root_segments[i], 1e-6f);
    }
    free(root_segments);
    return 1;
}

int sturm_solve4(real a, real b, real c, real d, real e, real *roots, uint *root_count)
{
    poly p;
    if (!poly_init(&p, 4))
        return 0;
    p.a[0] = a;
    p.a[1] = b;
    p.a[2] = c;
    p.a[3] = d;
    p.a[4] = e;
    poly_solve(&p, roots, root_count);
    poly_free(&p);
    return 1;
}

int real_root_count4(real a, real b, real c, real d, real e, uint *root_count)
{
    poly p;
    segment neg, pos;
    uint neg_count, pos_count;
    if (!poly_init(&p, 4))
        return 0;
    p.a[0] = a;
    p.a[1] = b;
    p.a[2] = c;
    p.a[3] = d;
    p.a[4] = e;
    poly_root_segments(&p, &neg, &pos);
    if (!poly_sturm_theorem(&p, &neg, &neg_count) ||
        !poly_sturm_theorem(&p, &pos, &pos_count)) {
        poly_free(&p);
        return 0;
    }
    poly_free(&p);
    *root_count = neg_count + pos_count;
    return 1;
}


















