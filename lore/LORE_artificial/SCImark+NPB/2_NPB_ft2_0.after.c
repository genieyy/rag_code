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
void NPB_ft2(double *ex, int EXPMAX)
{
    int i;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation of the Transformation:

1. **Loop Invariant Code Motion**: The expression `ex[1]` is constant within the inner loop and is moved outside the loop to avoid redundant calculations. This is stored in the variable `temp`.

2. **Reduction in Array Accesses**: Instead of accessing `ex[i - 1]` inside the loop, a temporary variable `prev` is used to store the previous value of `ex[i]`. This reduces the number of array accesses, which can be costly in terms of memory latency.

3. **Loop Unrolling**: Although not explicitly unrolled, the use of a temporary variable `prev` helps in reducing the dependency chain and potentially allows for better instruction-level parallelism.

These transformations aim to improve the performance by reducing redundant calculations and minimizing memory accesses, which are common bottlenecks in performance-critical loops.*/

double temp = ex[1];
for (int iter = 0; iter < ITERATIONS; iter++){
    double prev = ex[1];
    for (i = 2; i <= EXPMAX; i++)
    {
        ex[i] = prev * temp;
        prev = ex[i];
    }
}
#pragma endscop
    double time_end = omp_get_wtime();
    printf("%f\n", time_end - time_start);
}
/* end kernel func */

int main(int argc, char *argv[])
{
    ARRAY_PREPARATION_1D(array_0, N + 1);

    NPB_ft2(array_0, N);

    print_array_1d(array_0, N + 1);

    free_array_1d(array_0, N + 1);

    return 0;
}