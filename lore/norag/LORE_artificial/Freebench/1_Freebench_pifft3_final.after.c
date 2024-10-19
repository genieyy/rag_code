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
/*### Explanation:
1. **Type Casting to `double`**: By casting `radix` and `carry` to `double`, we avoid repeated casting within the loop, which can improve performance.
2. **Reduced Casting**: The `x` calculation is done using `double` to avoid multiple type conversions within the loop.
3. **Conditional Assignment**: The `carry` update is done using `double` to avoid repeated casting.
4. **Final Casting**: The final result is cast back to `int` to store in `out[j]`.
5. **Carry Update**: The `carry` variable is updated outside the loop to avoid repeated casting.

This transformation reduces the number of type conversions and improves the performance of the loop.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    double radix_d = (double)radix;
    double carry_d = (double)carry;
    for (int j = n - 1; j > 0; j--) {
        double x_d = (double)in1[j - 1] + (double)in2[j - 1] - carry_d;
        carry_d = (x_d >= radix_d) ? -1.0 : 0.0;
        out[j] = (int)(x_d - (radix_d * (int)carry_d));
    }
    carry = (int)carry_d;
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