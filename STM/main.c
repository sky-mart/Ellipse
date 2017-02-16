#include "ellipse_fitting.h"
#include "ellipse_fitting_tests.h" 

int main()
{
    int i;
       
    if (!testset_solve()) {
        i = 21;
    }
    if (!testset_dist_to_ellipse()) {
        i = 21;
    }
    if (!testset_linsolve()) {
        i = 21;
    }
    
    while (1) {}
}
