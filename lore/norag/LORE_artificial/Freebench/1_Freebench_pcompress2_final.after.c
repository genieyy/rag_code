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
/*### Explanation:
- **Reduced Division Operations**: By introducing a temporary variable `temp`, we reduce the number of division operations from two to one per iteration. This can provide a slight performance improvement, especially if division is more expensive than other operations.
- **Reuse of `temp`**: The temporary variable `temp` is reused to store the result of `(freq[i] + 1) / 2`, which is then used both to update `freq[i]` and to update `cum`. This avoids recalculating the same value multiple times.
- **No Redefinition of Variables**: The original variables are not redefined, and the new variable `temp` is of type `double`, as requested.

This version is expected to be slightly more efficient due to the reduced number of division operations and the reuse of the temporary variable.*/

for (int iter = 0; iter < ITERATIONS; iter++){
    double cum = 0.0;
    for (i = No_of_symbols; i > 0; i--) {
        double temp = (freq[i] + 1) / 2;
        freq[i] = temp;
        cum_freq[i] = cum;
        cum += temp;
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
