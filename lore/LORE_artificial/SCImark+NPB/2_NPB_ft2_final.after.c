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

3. **Loop Unrolling**: The inner loop is partially unrolled to reduce the number of iterations and improve instruction-level parallelism. This is done by multiplying the result by `temp` multiple times in a single iteration.

4. **Reduction in Loop Control**: By unrolling the loop, the number of loop control operations (incrementing `i` and checking the condition) is reduced, which can further improve performance.

These transformations aim to improve the performance by reducing redundant calculations, minimizing memory accesses, and reducing loop control overhead, which are common bottlenecks in performance-critical loops.*/

/*### Explanation of the Transformation:

1. **Loop Invariant Code Motion**: The expression `ex[1]` is constant within the inner loop and is moved outside the loop to avoid redundant calculations. This is stored in the variable `temp`.

2. **Reduction in Array Accesses**: Instead of accessing `ex[i - 1]` inside the loop, a temporary variable `prev` is used to store the previous value of `ex[i]`. This reduces the number of array accesses, which can be costly in terms of memory latency.

3. **Loop Unrolling**: The inner loop is partially unrolled to reduce the number of iterations and improve instruction-level parallelism. This is done by multiplying the result by `temp` multiple times in a single iteration.

4. **Reduction in Loop Control**: By unrolling the loop, the number of loop control operations (incrementing `i` and checking the condition) is reduced, which can further improve performance.

These transformations aim to improve the performance by reducing redundant calculations, minimizing memory accesses, and reducing loop control overhead, which are common bottlenecks in performance-critical loops.*/

double temp = ex[1];
for (int iter = 0; iter < ITERATIONS; iter++) {
    double prev = ex[1];
    for (i = 2; i <= EXPMAX - 3; i += 4) {
        ex[i] = prev * temp;
        prev = ex[i];
        ex[i + 1] = prev * temp;
        prev = ex[i + 1];
        ex[i + 2] = prev * temp;
        prev = ex[i + 2];
        ex[i + 3] = prev * temp;
        prev = ex[i + 3];
    }
    for (; i <= EXPMAX; i++) {
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