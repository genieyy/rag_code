#include <stdio.h>  
#include <stdlib.h>  
#include <string.h>
#include <math.h>
#include "lore.h"

#include <omp.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y) ((x) > (y)? (x) : (y))
#define min(x,y) ((x) < (y)? (x) : (y))

/* start param define */
#define N 20000
#define M 13
/* end parameters define */

/* start kernel func */
void NPB_bt9(double *dtemp, double ** ce, int n_, int m_) {
    int i, j, k, m;
    double dssp = 0.1; // Example value, replace with actual value  
    double xi = 5, eta = 14, zeta = 16; // Assuming these are defined elsewhere  

    double time_start = omp_get_wtime();
#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int m = 0; m < n_; m += 16) {
        for (int i = m; i < m + 16 && i < n_; i++) {
            dtemp[i] = ce[i][0] +
                xi * (ce[i][1] + xi * (ce[i][4] + xi * (ce[i][7]
                    + xi * ce[i][10]))) +
                eta * (ce[i][2] + eta * (ce[i][5] + eta * (ce[i][8]
                    + eta * ce[i][11]))) +
                zeta * (ce[i][3] + zeta * (ce[i][6] + zeta * (ce[i][9] +
                    zeta * ce[i][12])));
        }
    }
}
#pragma endscop
    double time_end = omp_get_wtime();
    printf("%f\n", time_end - time_start);
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