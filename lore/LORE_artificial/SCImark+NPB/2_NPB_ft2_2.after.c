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
/*### Explanation of Optimizations:
1. **Loop Invariant Code Motion**: The calculation of `ex[1]` is moved outside the inner loop because it does not change within the loop. This reduces redundant calculations.
2. **Reduction in Array Accesses**: Instead of accessing `ex[i - 1]` repeatedly, a temporary variable `current` is used to store the value of the previous iteration's result. This reduces the number of array accesses, which can be costly.
3. **Reduction in Multiplications**: By using the `current` variable, the number of multiplications is reduced to one per iteration of the inner loop, instead of two (one for `ex[i - 1]` and one for `ex[1]`).

These optimizations help improve the performance of the loop by reducing redundant calculations and minimizing the number of array accesses.*/

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