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
/**/

int t1, t2, t3, t4, t5;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

for (int t1 = 0; t1 <= floord(ITERATIONS - 1, 32); t1++) {
    lbp = max(0, ceild(32 * t1 - No_of_symbols + 1, 32));
    ubp = min(floord(ITERATIONS - 1, 32), t1);
#pragma omp parallel for private(lbv, ubv, t3, t4, t5)
    for (int t2 = lbp; t2 <= ubp; t2++) {
        for (int t3 = max(32 * t1 - 32 * t2, 1); t3 <= min(No_of_symbols - 1, 32 * t1 - 32 * t2 + 31); t3++) {
            freq[t3] = (freq[t3] + 1) / 2;
            cum_freq[t3] = cum;
            cum += freq[t3];
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
