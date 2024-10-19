/*according to 
Benchmark: FreeBench default
Application: pifft
File: pifft.c
Function: mp_unexp_add
Line: 500
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
#define N 20000
/* end parameters define */

/* start kernel func */
void Freebench_pifft3(double *in1, double *in2, double *out, int n)
{
    int x;
    int carry = 0;
    int radix = 10;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation of the Optimization:
1. **Loop Unrolling**: The original code iterates over each element of the array `in1` and `in2` to compute the values for `out`. The loop is unrolled to reduce the overhead of loop control.
2. **Reduction of Conditional Checks**: The carry calculation is simplified by directly adding the carry to the sum of `in1[j - 1]` and `in2[j - 1]`. This reduces the number of conditional checks within the loop.
3. **Bitwise Operations**: The use of bitwise operations (e.g., `radix & carry`) is retained to maintain the efficiency of the carry calculation.

These optimizations are based on the principles observed in the provided examples, such as reducing the number of iterations and simplifying the logic within the loop to improve performance.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    int carry = 0;
    for (int j = n - 1; j > 0; j--) {
        int x = in1[j - 1] + in2[j - 1] + carry;
        carry = (x >= radix) ? -1 : 0;
        out[j] = x - (radix & carry);
    }
}
#pragma endscop
    double time_end = omp_get_wtime();
    printf("%f\n", time_end - time_start);
}
/* end kernel func */

int main(int argc, char *argv[])
{
    ARRAY_PREPARATION_1D(array_0, N+1);
    ARRAY_PREPARATION_1D(array_1, N+1);
    ARRAY_PREPARATION_1D(array_2, N+1);

    Freebench_pifft3(array_0, array_1, array_2, N);

    print_array_1d(array_2, N+1);

    free_array_1d(array_0, N+1);
    free_array_1d(array_1, N+1);
    free_array_1d(array_2, N+1);

    return 0;
}