//
//  main.c
//  Eq4SinglePrec
//
//  Created by Влад Агиевич on 27.02.17.
//  Copyright © 2017 Mocsmart. All rights reserved.
//

#include <stdio.h>
#include "equations.h"

int main(int argc, const char * argv[])
{
    real a, b, x[2], d, l;
    
    a = 1.50000203f;
    b = 0.999994993f;
    x[0] = 0.655670345f;
    x[1] = 0.899404049f;
    
    dist_to_ellipse(a, b, x, &d, &l);
    
    printf("d = %e, l = %e\n", d, l);
    
    return 0;
}
