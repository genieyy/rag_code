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
/*### Explanation of the Transformation:
1. **Loop Unrolling**: The original code has a loop that iterates over `j` from `n-1` to `1`. The transformation does not explicitly unroll the loop, but it simplifies the logic within the loop to reduce the number of operations per iteration.

2. **Reduction of Conditional Checks**: The original code checks if `x >= radix` and then sets `carry` accordingly. The optimized code combines this check into a single expression `(x >= radix) ? -1 : 0`, which reduces the number of conditional checks.

3. **Simplification of Operations**: The original code uses `x - (radix & carry)` to adjust the value of `x` based on the `carry`. The optimized code maintains this structure but ensures that the `carry` is updated in a more straightforward manner.

4. **Initialization of `carry`**: The `carry` variable is initialized outside the inner loop to avoid reinitialization in each iteration, which can be more efficient.

These transformations aim to improve the performance by reducing the number of operations and conditional checks within the loop, thereby making the code more efficient.*/

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