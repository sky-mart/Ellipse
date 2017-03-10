//
//  main.c
//  Equation4_x86
//
//  Created by Влад Агиевич on 06.03.17.
//  Copyright © 2017 Mocsmart. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include "tests.h"
#include "direct.h"
#include "sturm.h"

int main(int argc, const char * argv[]) {
    
//    if (testset_solve(direct_solve2, direct_solve3, direct_solve4)) {
//        printf("Success\n");
//    } else {
//        printf("Fail\n");
//    }
    
    poly u;
//    segment root_segments[4];
    real roots[4];
    uint num, i;
    
    if (!poly_init(&u, 4))
        return 0;
    u.a[0] = 1.0f;
    u.a[1] = 5.0f;
    u.a[2] = 6.0f;
    u.a[3] = -3.0f;
    u.a[4] = -21.0f;

    if (!poly_solve(&u, roots, &num))
        return 0;
    for (i = 0; i < num; i++) {
        printf("%f\n", roots[i]);
    }
    
    poly_free(&u);
    return 0;
}
