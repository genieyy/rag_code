/*according to 
Benchmark: FreeBench default
Application: pifft
File: pifft.c
Function: mp_sprintf
Line: 1450
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
#define N 200
#define M 300
/* end parameters define */

/* start kernel func */
void Freebench_pifft2(double *in, double *out, int n, int log10_radix)
{
    int x, y;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation:
1. **Loop Tiling and Parallelization**: The outer loop is tiled and parallelized using OpenMP to maximize the use of multiple cores.
2. **Reduced Array Accesses**: The `offset` variable is precomputed to reduce the number of arithmetic operations inside the innermost loop.
3. **Register Usage**: The `lbv` and `ubv` variables are marked as `register` to suggest the compiler to use CPU registers for these variables, potentially improving performance.
4. **Bounds Checking**: The inner loop bounds are adjusted to ensure that the loop only processes valid indices, reducing unnecessary iterations.

This version should provide a further optimized performance compared to the previous versions.*/

int t1, t2, t3, t4;
int lb, ub, lbp, ubp;
register int lbv, ubv;

for (int t1 = 0; t1 <= floord(ITERATIONS, 32); t1++) {
    lbp = max(0, t1 - 1);
    ubp = min(floord(ITERATIONS, 32), t1);
#pragma omp parallel for private(lbv, ubv, t3, t4)
    for (int t2 = lbp; t2 <= ubp; t2++) {
        for (int t3 = max(3, 32 * t1); t3 <= min(n + 1, 32 * t1 + 31); t3++) {
            x = in[t3];
            int offset = (t3 - 3) * 3;
            for (int t4 = log10_radix - 1; t4 >= 0; t4--) {
                y = x % 10;
                x /= 10;
                out[offset + t4] = (48 + y);
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
    ARRAY_PREPARATION_1D(array_0, N+2);

    ARRAY_PREPARATION_1D(array_1, 3*N+M+1);

    Freebench_pifft2(array_0, array_1, N, M);

    print_array_1d(array_1, 3*N+M+1);

    free_array_1d(array_0, N+2);
    
    free_array_1d(array_1, 3*N+M+1);

    return 0;
}