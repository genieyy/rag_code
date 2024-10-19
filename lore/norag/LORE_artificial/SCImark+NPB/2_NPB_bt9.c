#include <stdio.h>  
#include <stdlib.h>  
#include <string.h>
#include <math.h>
#include "lore.h"

/* start param define */
#define N 20000
#define M 13
/* end parameters define */

/* start kernel func */
void NPB_bt9(double *dtemp, double ** ce, int n_, int m_) {
    int i, j, k, m;
    double dssp = 0.1; // Example value, replace with actual value  
    double xi = 5, eta = 14, zeta = 16; // Assuming these are defined elsewhere  

#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++){
    // Your provided code here  
    for (m = 0; m < n_; m++) {
        dtemp[m] = ce[m][0] +
            xi * (ce[m][1] + xi * (ce[m][4] + xi * (ce[m][7]
                + xi * ce[m][10]))) +
            eta * (ce[m][2] + eta * (ce[m][5] + eta * (ce[m][8]
                + eta * ce[m][11]))) +
            zeta * (ce[m][3] + zeta * (ce[m][6] + zeta * (ce[m][9] +
                zeta * ce[m][12])));
    }
}
#pragma endscop
}
/* end kernel func */

int main() {
    
    ARRAY_PREPARATION_1D(array_0, N+1);
    ARRAY_PREPARATION_2D(array_1, N+1, M+1);

    NPB_bt9(array_0, array_1, N, M);
 
    print_array_1d(array_0, N+1);

    free_array_1d(array_0, N+1);
    free_array_2d(array_1, N+1, M+1);
    return 0;
}