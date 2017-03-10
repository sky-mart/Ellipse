//
//  direct.h
//  Equation4_x86
//
//  Direct solution of 4th degree equation
//
//  Created by Влад Агиевич on 06.03.17.
//  Copyright © 2017 Mocsmart. All rights reserved.
//

#ifndef direct_h
#define direct_h

#include "settings.h"

int direct_solve2(dreal a, dreal b, dreal c, dcomplex roots[2]);
int direct_solve3(dreal a, dreal b, dreal c, dreal d, dcomplex roots[3]);
int direct_solve4(dreal a, dreal b, dreal c, dreal d, dreal e, dcomplex roots[4]);

int direct_solve4_real(dreal a, dreal b, dreal c, dreal d, dreal e, real *roots, uint *num_roots);

#endif /* direct_h */
