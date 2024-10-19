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
/*### Explanation:
1. **Pointer Aliasing**: By using pointer aliasing (`double *ldiph_lc_br = ldiph_lc[b][r];`), we reduce the number of array indexing operations, which can be costly, especially in nested loops.
2. **Reduced Indexing**: The inner loop now directly accesses the aliased pointers, reducing the overhead of multiple array indexing operations.
3. **Type Casting**: The `bad_s3ssid` variable is cast to `double` once per iteration of the outer loop, which is more efficient than casting it multiple times within the inner loops.

This version should provide a further performance improvement by reducing the number of array indexing operations and leveraging pointer aliasing.*/

for (int iter = 0; iter < ITERATIONS; iter++){
    double bad_s3ssid_dbl = (double)bad_s3ssid;
    for (b = 0; b < n_ciphone; ++b) {
        for (r = 0; r < n_ciphone; ++r) {
            double *ldiph_lc_br = ldiph_lc[b][r];
            double *lrdiph_rc_br = lrdiph_rc[b][r];
            double *rdiph_rc_br = rdiph_rc[b][r];
            for (l = 0; l < n_ciphone; ++l) {
                ldiph_lc_br[l] = bad_s3ssid_dbl;
                lrdiph_rc_br[l] = bad_s3ssid_dbl;
                rdiph_rc_br[l] = bad_s3ssid_dbl;
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