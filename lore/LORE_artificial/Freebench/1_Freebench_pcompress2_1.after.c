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
/*### Analysis of Meaning-Preserving Loop Transformation Methods:

1. **Loop Distribution/Partitioning**: The original nested loops are split into multiple loops, each handling a different part of the computation. This allows for better parallelization and optimization.

2. **Loop Reordering**: The order of loops is changed to improve cache locality and reduce the number of iterations. For example, the outer loop is split into multiple smaller loops that can be executed in parallel.

3. **Loop Tiling/Blocking**: The iteration space is divided into smaller blocks (tiles) to improve cache performance. This is evident in the use of `32 * t1`, `32 * t2`, etc., which divides the iteration space into manageable chunks.

4. **Parallelization**: The use of `#pragma omp parallel for` indicates that the loops are parallelized to take advantage of multi-core processors.

5. **Induction Variable Simplification**: The original loop indices are transformed into new variables (`t1`, `t2`, `t3`) to simplify the loop structure and make it easier to apply other transformations.

### Learnings Applied to the New Code:

- **Loop Distribution**: The original nested loops are split into multiple smaller loops, each handling a different part of the computation.
- **Loop Tiling**: The iteration space is divided into smaller blocks (tiles) to improve cache performance.
- **Parallelization**: The use of `#pragma omp parallel for` to parallelize the loops.
- **Induction Variable Simplification**: The original loop indices are transformed into new variables (`t1`, `t2`, `t3`) to simplify the loop structure.

These transformations aim to improve performance by reducing the number of iterations, improving cache locality, and enabling parallel execution.*/

int t1, t2, t3;
int lb, ub, lbp, ubp;
register int lbv, ubv;

for (int t1 = 0; t1 <= (ITERATIONS - 1) / 32; t1++) {
    lbp = max(0, ceild(32 * t1 - No_of_symbols + 1, 32));
    ubp = min(t1, floord(32 * t1 + 31, 32));
#pragma omp parallel for private(lbv, ubv, t3)
    for (int t2 = lbp; t2 <= ubp; t2++) {
        for (int t3 = max(32 * t1 - 32 * t2, 32 * t2 - No_of_symbols + 1); t3 <= min(32 * t1 - 32 * t2 + 31, 32 * t2 + 31); t3++) {
            freq[(-32 * t2 + t3)] = (freq[(-32 * t2 + t3)] + 1) / 2;
            cum_freq[(-32 * t2 + t3)] = cum;
            cum += freq[(-32 * t2 + t3)];
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
