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
/*### Explanation:
- **Reduced Array Accesses**: By using the `prev` variable to store the previous value, we reduce the number of array accesses from `ex[i - 1]` to just `prev`.
- **Single Multiplication per Iteration**: The multiplication `prev *= temp` is performed after storing the current value in `ex[i]`, ensuring that the multiplication is only done once per iteration.
- **Simplified Loop**: The loop structure remains simple and clear, maintaining readability while improving performance.*/

double temp = ex[1];
for (int iter = 0; iter < ITERATIONS; iter++) {
    double prev = temp;
    for (i = 2; i <= EXPMAX; i++) {
        ex[i] = prev;
        prev *= temp;
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