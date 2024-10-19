#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lore.h"
#include <omp.h>
/*
Cortexsuite\cortex\sphinx\pocketsphinx\dict2pid.c
line 422 - line 430
*/

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
  int t1, t2, t3, t4, t5, t6, t7, t8;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
lbp=0;
ubp=floord(n_ciphone-1,32);
#pragma omp parallel for private(lbv,ubv,t2,t3,t4,t5,t6,t7,t8)
for (t1=lbp;t1<=ubp;t1++) {
  for (t2=0;t2<=floord(n_ciphone-1,32);t2++) {
    for (t3=0;t3<=floord(n_ciphone-1,32);t3++) {
      for (t4=0;t4<=floord(ITERATIONS-1,32);t4++) {
        for (t5=32*t1;t5<=min(n_ciphone-1,32*t1+31);t5++) {
          for (t6=32*t2;t6<=min(n_ciphone-1,32*t2+31);t6++) {
            for (t7=32*t3;t7<=min(n_ciphone-1,32*t3+31);t7++) {
              for (t8=32*t4;t8<=min(ITERATIONS-1,32*t4+31);t8++) {
                ldiph_lc[t5][t6][t7] = bad_s3ssid;;
                lrdiph_rc[t5][t7][t6] = bad_s3ssid;;
                rdiph_rc[t5][t7][t6] = bad_s3ssid;;
              }
            }
          }
        }
      }
    }
  }
}
/* End of CLooG code */
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
