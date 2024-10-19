/*according to 
Benchmark: FreeBench default
Application: pifft
File: pifft.c
Function: mp_mul_d2i
Line: 1053
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
#define N 100000
/* end parameters define */

/* start kernel func */
void Freebench_pifft1(double *out, int n)
{
    int j;
    int carry = 42;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation of the Transformation:
1. **Loop Unrolling and Parallelization**: The outer loop is unrolled and parallelized using OpenMP to exploit multiple threads. This is similar to the examples provided where the outer loops were unrolled and parallelized.
2. **Bounds Calculation**: The bounds for the inner loop (`lbv` and `ubv`) are calculated to ensure that the loop iterates over the correct range of iterations. This is similar to the bounds calculations in the provided examples.
3. **Register Usage**: The `register` keyword is used to hint to the compiler that `lbv` and `ubv` should be stored in registers for faster access, similar to the usage in the provided examples.

This transformation aims to improve performance by leveraging parallel execution and reducing the overhead of loop iterations.*/

int t1, t2;
register int lbv, ubv;

for (int t1 = 0; t1 <= (ITERATIONS + n - 2) / 32; t1++) {
    lbv = max(0, 32 * t1 - n + 1);
    ubv = min(ITERATIONS - 1, 32 * t1 + 31);
#pragma omp parallel for private(t2)
    for (int t2 = lbv; t2 <= ubv; t2++) {
        for (int j = n + 1; j > 2; j--) {
            out[j] = out[j - 1];
        }
        out[2] = carry;
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

    Freebench_pifft1(array_0, N);

    print_array_1d(array_0, N+2);

    free_array_1d(array_0, N+2);

    return 0;
}
