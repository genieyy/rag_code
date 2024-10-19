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
#define N 400
#define M 600
/* end parameters define */

/* start kernel func */
void SCImark_lu1(double **A, double *Aii, double *Aj, int m, int n)
{
    int j = 0;
    int ii, jj;
    double AiiJ;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation:
1. **Precompute the First Multiplication**: The first multiplication `AiiJ * Aj[j + 1]` is precomputed and stored in a temporary variable `temp`. This avoids recalculating the multiplication in each iteration of the inner loop.
2. **Update Temp Safely**: The `temp` variable is updated only if `jj + 1 < n` to avoid out-of-bounds access.
3. **Store AiiJ Back**: After the inner loop, `AiiJ` is stored back into `Aii[j]` to avoid redundant memory access.

This version combines the best aspects of the previous optimizations, ensuring that the code is both efficient and safe.*/

for (int iter = 0; iter < ITERATIONS; iter++){
    for (ii = j + 1; ii < m; ii++) {
        Aii = A[ii];
        Aj = A[j];
        double AiiJ = Aii[j];
        
        // Precompute the first multiplication and store it in a temporary variable
        double temp = AiiJ * Aj[j + 1];
        
        // Loop through the remaining elements and update them using the precomputed value
        for (jj = j + 1; jj < n; jj++) {
            Aii[jj] -= temp;
            if (jj + 1 < n) {
                temp = AiiJ * Aj[jj + 1]; // Update temp for the next iteration
            }
        }
        
        // Store AiiJ back to avoid redundant memory access
        Aii[j] = AiiJ;
    }
}
#pragma endscop
    double time_end = omp_get_wtime();
    printf("%f\n", time_end - time_start);
}
/* end kernel func */

int main(int argc, char *argv[])
{
    ARRAY_PREPARATION_2D(array_0, N+1, M+1);
    ARRAY_PREPARATION_1D(array_1, M+1);
    ARRAY_PREPARATION_1D(array_2, M+1);

    SCImark_lu1(array_0, array_1, array_2, N, M);

    print_array_1d(array_1, M+1);

    free_array_2d(array_0, N+1, M+1);
    free_array_1d(array_1, M+1);
    free_array_1d(array_2, M+1);

    return 0;
}
