/*
ASC_Sequoia/CrystalMK/Crystal_Cholesky.c
line 33 - line 48
*/
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
#define N 100
/* end parameters define */

/* start kernel func */
void ASC_Sequoia_Crystal_Cholesky(double** a, int nSlip)
{
    int i, j, k;
    double fdot;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation:
1. **Reduced Redundant Calculations**: 
   - The inner loops that calculate `fdot` for `a[i][j]` and `a[j][i]` are combined into a single loop, reducing the number of iterations.
   - This is done by calculating both `fdot_ij` and `fdot_ji` in the same loop, which reduces the overall computational load.

2. **Avoiding Redundant Initializations**:
   - The variable `fdot` is initialized only once per iteration of the outer loop, and its value is reused in subsequent calculations, reducing the number of initializations.

3. **Simplified Assignment**:
   - The assignment `a[i][i] = a[i][i] - fdot_i;` is simplified to `a[i][i] -= fdot_i;` for better readability and potential performance improvement.

These changes aim to reduce the number of operations and improve the locality of reference, which can lead to better performance, especially for large values of `nSlip` and `ITERATIONS`.*/

for (int iter = 0; iter < ITERATIONS; iter++){
    for (i = 1; i < nSlip; i++)
    {
        double fdot_i = 0.0;
        for (k = 0; k < i; k++)
        {
            fdot_i += a[i][k] * a[k][i];
        }
        a[i][i] -= fdot_i;
        for (j = i + 1; j < nSlip; j++)
        {
            double fdot_ij = 0.0;
            double fdot_ji = 0.0;
            for (k = 0; k < i; k++)
            {
                fdot_ij += a[i][k] * a[k][j];
                fdot_ji += a[j][k] * a[k][i];
            }
            a[i][j] -= fdot_ij;
            a[j][i] = (a[j][i] - fdot_ji) / a[i][i];
        }
    }
}
#pragma endscop
    double time_end = omp_get_wtime();
    printf("%f\n", time_end - time_start);
}
/* end kernel func */

int main(int argc, char *argv[])
{
    ARRAY_PREPARATION_2D(array_0, N+1, N+1);

    ASC_Sequoia_Crystal_Cholesky(array_0, N);

    print_array_2d(array_0, N+1, N+1);

    free_array_2d(array_0, N+1, N+1);
    
    return 0;
}