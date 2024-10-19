/*according to 
Benchmark: FreeBench default
Application: pifft
File: pifft.c
Function: mp_mul_d2i
Line: 1053
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
#define N 100000
/* end parameters define */

/* start kernel func */
void Freebench_pifft1(double *out, int n)
{
    int j;
    int carry = 42;

    double time_start = omp_get_wtime();
#pragma scop
/**/

int t1, t2, t3;
int lb, ub, lbp, ubp;
register int lbv, ubv;

for (int t1 = 0; t1 <= floord(ITERATIONS + n - 2, 32); t1++) {
    lbp = max(0, ceild(32 * t1 - n + 1, 32));
    ubp = min(floord(ITERATIONS - 1, 32), t1);
#pragma omp parallel for private(lbv, ubv, t3)
    for (int t2 = lbp; t2 <= ubp; t2++) {
        for (int t3 = max(32 * t1 - 32 * t2, 3); t3 <= min(n + 1, 32 * t1 - 32 * t2 + 31); t3++) {
            lbv = max(32 * t2, t3 - 1);
            ubv = min(32 * t2 + 31, ITERATIONS - 1);
#pragma ivdep
#pragma vector always
            for (int iter = lbv; iter <= ubv; iter++) {
                out[t3] = out[t3 - 1];
            }
        }
        if (t1 == t2) {
            lbv = max(32 * t2, 2);
            ubv = min(32 * t2 + 31, ITERATIONS - 1);
#pragma ivdep
#pragma vector always
            for (int iter = lbv; iter <= ubv; iter++) {
                out[2] = carry;
            }
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
    ARRAY_PREPARATION_1D(array_0, N+2);

    Freebench_pifft1(array_0, N);

    print_array_1d(array_0, N+2);

    free_array_1d(array_0, N+2);

    return 0;
}
