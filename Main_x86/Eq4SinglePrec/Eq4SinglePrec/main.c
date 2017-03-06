//
//  main.c
//  Eq4SinglePrec
//
//  Created by Влад Агиевич on 27.02.17.
//  Copyright © 2017 Mocsmart. All rights reserved.
//

#include <stdio.h>
//#include "ellipse_fitting.h"
//
//typedef union {
//    float f;
//    struct {
//        unsigned int mantisa : 23;
//        unsigned int exponent : 8;
//        unsigned int sign : 1;
//    } parts;
//} double_cast;
//
//void float_repr(float x)
//{
//    uint i = (*(uint *)(&x));
//    uint8_t s = i >> 31;
//    int8_t e = ((i & 0x7f800000) >> 23) - 127;
//    int m = (i & 0x007fffff) | 0x00800000;
//    printf("i=%x\n", i);
//    printf("s=%i\n", s);
//    printf("e=%i\n", e);
//    printf("m=%x\n", m);
//    printf("repr=%f\n", ((float)m) * (1 << e) / (1 << 23));
//}
//
//void print_float_hex(float x)
//{
//    uint i = (*(uint *)(&x));
//    printf("%x\n", i);
//}
//
//void print_float_bin(float x)
//{
//    uint itg = (*(uint *)(&x));
//    
//    int i;
//    for (i = 31; i >= 0; i--) {
//        printf("%i", (itg & (1 << i)) >> i);
//        if ((i == 31) || (i == 23))
//            printf(" ");
//    }
//    printf("\n");
//}
//
//int main(int argc, const char * argv[])
//{
//    real a, b, x[2], d, l;
//    
//    a = 1.50000203f;
//    b = 0.999994993f;
//    x[0] = 0.655670345f;
//    x[1] = 0.899404049f;
//    
//    dist_to_ellipse(a, b, x, &d, &l);
//
//    printf("d = %e, l = %e\n", d, l);
////    double_cast d1;
////    d1.f = 4.2f;
////    printf("sign = %i\n",d1.parts.sign);
////    printf("exponent = %i\n",d1.parts.exponent);
////    printf("mantisa = %i\n",d1.parts.mantisa);
////    printf("%i\n", 1 << 2);
////    float_repr(3.2f);
////    int i = 0;
////    for (i = 0; i < 10; i++) {
////        print_float_bin((float)(1 << i));
////    }
//
//    return 0;
//}

#include "ellipse.h"
#include "ellipse_fitting_tests.h"

extern uint solve3_call_counter;
extern uint single_prec_not_enough_counter;

int main()
{
    int i;
    
    //srand(time(NULL));
    
    if (!testset_solve()) {
        i = 21;
        printf("testset_solve: failed\n");
    } else {
        printf("testset_solve: passed\n");
    }
    if (!testset_dist_to_ellipse()) {
        i = 21;
        printf("testset_dist_to_ellipse: failed\n");
    } else {
        printf("testset_dist_to_ellipse: passed\n");
    }
    if (!testset_linsolve()) {
        i = 21;
        printf("testset_linsolve: failed\n");
    } else {
        printf("testset_linsolve: passed\n");
    }
    if (!testset_ellipse_fitting()) {
        i = 21;
        printf("testset_ellipse_fitting: failed\n");
    } else {
        printf("testset_ellipse_fitting: passed\n");
    }
    
    printf("solve3_call_counter: %i\n", solve3_call_counter);
    printf("single_prec_not_enough_counter: %i\n", single_prec_not_enough_counter);
    
//    while (1) {}
    return 1;
}

