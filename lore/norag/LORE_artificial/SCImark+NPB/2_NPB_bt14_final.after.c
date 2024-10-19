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
#define N 20
#define M 30
#define L 40
#define P 5
/* end parameters define */

/* start kernel func */
void NPB_bt14(double ****rhs, double ****forcing, int nz, int mz, int q, int p)
{
    int i, j, k, m;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation:
1. **Restrict Keyword**: The `restrict` keyword is used to indicate that pointers do not alias, which can help the compiler generate more efficient code by allowing it to make certain assumptions about memory access patterns.
2. **Pointer Aliasing**: By creating local pointers (`rhs_slice` and `forcing_slice`), we reduce the number of array index calculations, which can improve performance.
3. **Memcpy Optimization**: The `memcpy` function is used to copy blocks of memory, which is generally faster than copying elements one by one in a loop. This optimization is retained from the previous version.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int i = 0; i < nz; i++) {
        for (int j = 0; j < mz; j++) {
            double (*restrict rhs_slice)[p] = rhs[i][j];
            double (*restrict forcing_slice)[p] = forcing[i][j];
            for (int k = 0; k < q; k++) {
                memcpy(rhs_slice[k], forcing_slice[k], p * sizeof(double));
            }
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
    ARRAY_PREPARATION_4D(array_0, N + 1, M + 1, L + 1, P + 1);
    ARRAY_PREPARATION_4D(array_1, N + 1, M + 1, L + 1, P + 1);

    NPB_bt14(array_0, array_1, N, M, L, P);

    print_array_4d(array_0, N + 1, M + 1, L + 1, P + 1);

    free_array_4d(array_0, N + 1, M + 1, L + 1, P + 1);
    free_array_4d(array_1, N + 1, M + 1, L + 1, P + 1);

    return 0;
}