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

lbp = 0;
ubp = floord(n_ciphone - 1, 32);

#pragma omp parallel for private(lbv, ubv, t2, t3, t4)
for (t1 = lbp; t1 <= ubp; t1++) {
    for (t2 = 0; t2 <= ubp; t2++) {
        for (t3 = 0; t3 <= ubp; t3++) {
            for (t4 = max(0, 32 * t1); t4 <= min(n_ciphone - 1, 32 * t1 + 31); t4++) {
                for (int b = max(0, 32 * t2); b <= min(n_ciphone - 1, 32 * t2 + 31); b++) {
                    for (int r = max(0, 32 * t3); r <= min(n_ciphone - 1, 32 * t3 + 31); r++) {
                        for (int l = max(0, 32 * t3); l <= min(n_ciphone - 1, 32 * t3 + 31); l++) {
                            ldiph_lc[b][r][l] = bad_s3ssid;
                            lrdiph_rc[b][l][r] = bad_s3ssid;
                            rdiph_rc[b][l][r] = bad_s3ssid;
                        }
                    }
                }
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