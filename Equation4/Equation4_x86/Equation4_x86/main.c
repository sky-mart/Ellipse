//
//  main.c
//  Equation4_x86
//
//  Created by Влад Агиевич on 06.03.17.
//  Copyright © 2017 Mocsmart. All rights reserved.
//

#include <stdio.h>
#include "tests.h"
#include "direct.h"

int main(int argc, const char * argv[]) {
    // insert code here...
    if (testset_solve(direct_solve2, direct_solve3, direct_solve4)) {
        printf("Success\n");
    } else {
        printf("Fail\n");
    }
    return 0;
}
