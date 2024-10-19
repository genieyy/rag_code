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

1. **Loop Invariant Code Motion**: The expression `ex[1]` is computed outside the inner loop because it does not change within the loop. This reduces the number of redundant computations.

2. **Strength Reduction**: Instead of repeatedly multiplying `ex[i - 1]` by `ex[1]`, we use a temporary variable `current` to store the intermediate result of the exponentiation. This reduces the number of array accesses and multiplications.

3. **Reduction in Array Accesses**: By using the `current` variable, we avoid accessing `ex[i - 1]` in each iteration of the inner loop, which can be more efficient, especially if `ex` is stored in memory that is not cache-friendly.

These transformations aim to improve performance by reducing redundant computations and minimizing memory accesses.*/

double temp = ex[1];
for (int iter = 0; iter < ITERATIONS; iter++){
    double current = temp;
    for (i = 2; i <= EXPMAX; i++)
    {
        ex[i] = current;
        current *= temp;
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