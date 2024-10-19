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
void Freebench_pcompress2(double *freq, double *cum_freq, int No_of_symbols)
{
    int i;
    int cum = 0;

    double time_start = omp_get_wtime();
#pragma scop
/*The examples provided demonstrate several loop transformation techniques that can be used to optimize performance. These techniques include loop tiling, loop unrolling, and parallelization using OpenMP. Here's a summary of the key methods used:

1. **Loop Tiling**: This technique breaks the iteration space into smaller blocks (tiles) to improve cache locality and reduce cache misses.
2. **Loop Unrolling**: This technique reduces loop overhead by performing multiple operations within a single iteration of the loop.
3. **Parallelization**: Using OpenMP to parallelize loops can significantly improve performance by distributing the workload across multiple threads.

Applying these techniques to the given code:



### Explanation:
1. **Parallelization**: The outer loop is parallelized using OpenMP to distribute the iterations across multiple threads.
2. **Loop Unrolling**: The inner loop is unrolled by a factor of 2 to reduce loop overhead and improve instruction-level parallelism.
3. **Loop Ordering**: The inner loop iterates from `No_of_symbols - 1` to `0` to maintain the original loop semantics while allowing for unrolling.

This optimized code should provide better performance by leveraging parallel execution and reducing loop overhead through unrolling.*/

#pragma omp parallel for private(i, cum)
for (int iter = 0; iter < ITERATIONS; iter++) {
    double cum = 0.0;
    for (int i = No_of_symbols - 1; i >= 0; i -= 2) { // Loop unrolling by 2
        freq[i] = (freq[i] + 1) / 2;
        cum_freq[i] = cum;
        cum += freq[i];

        if (i > 0) { // Unroll the next iteration
            freq[i - 1] = (freq[i - 1] + 1) / 2;
            cum_freq[i - 1] = cum;
            cum += freq[i - 1];
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
    ARRAY_PREPARATION_1D(array_0, N+1);
    ARRAY_PREPARATION_1D(array_1, N+1);

    Freebench_pcompress2(array_0, array_1, N);

    print_array_1d(array_0, N+1);

    free_array_1d(array_0, N+1);

    return 0;
}
