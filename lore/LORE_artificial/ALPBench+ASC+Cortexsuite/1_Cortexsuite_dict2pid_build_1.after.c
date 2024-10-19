/*
Cortexsuite\cortex\sphinx\pocketsphinx\dict2pid.c
line 422 - line 430
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
#define N 40
/* end parameters define */

/* start kernel func */
void Cortexsuite_dict2pid_build(double ***ldiph_lc, double ***lrdiph_rc, double ***rdiph_rc, int n_ciphone)
{
    int b, r, l;
    double bad_s3ssid = 1;

    double time_start = omp_get_wtime();
#pragma scop
/**/

int t1, t2, t3, t4;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

for (int t1 = 0; t1 <= floord(ITERATIONS - 1, 16); t1++) {
    lbp = max(0, ceild(32 * t1 - ITERATIONS + 1, 32));
    ubp = floord(t1, 2);
#pragma omp parallel for private(lbv, ubv, t3, t4)
    for (int t2 = lbp; t2 <= ubp; t2++) {
        for (int t3 = max(32 * t1 - 32 * t2, 32 * t2 + 1); t3 <= min(n_ciphone - 1, 32 * t1 - 32 * t2 + 31); t3++) {
            lbv = 32 * t2;
            ubv = min(32 * t2 + 31, t3 - 1);
#pragma ivdep
#pragma vector always
            for (int t4 = lbv; t4 <= ubv; t4++) {
                ldiph_lc[t3][t4][t4] = bad_s3ssid;
                lrdiph_rc[t3][t4][t4] = bad_s3ssid;
                rdiph_rc[t3][t4][t4] = bad_s3ssid;
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
    ARRAY_PREPARATION_3D(array_0, N+1, N+1, N+1);
    ARRAY_PREPARATION_3D(array_1, N+1, N+1, N+1);
    ARRAY_PREPARATION_3D(array_2, N+1, N+1, N+1);

    Cortexsuite_dict2pid_build(array_0, array_1, array_2, N);

    print_array_3d(array_2, N+1, N+1, N+1);

    free_array_3d(array_0, N+1, N+1, N+1);
    free_array_3d(array_1, N+1, N+1, N+1);
    free_array_3d(array_2, N+1, N+1, N+1);
    
    return 0;
}