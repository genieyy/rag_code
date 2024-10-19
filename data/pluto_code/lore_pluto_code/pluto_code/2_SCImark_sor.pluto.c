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
#define M 50
#define L 60
/* end parameters define */

/* start kernel func */
void SCImark_sor(double **G, double *Gi, double *Gim1, double *Gip1, int num_iterations, int Mm1, int Nm1)
{
    int p, i, j;
	double one_minus_omega = 12;
    double omega_over_four = 48;

    double time_start = omp_get_wtime();
  int t1, t2, t3, t4, t5, t6;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
for (t1=0;t1<=ITERATIONS-1;t1++) {
  for (t2=0;t2<=num_iterations-1;t2++) {
    for (t3=0;t3<=floord(Mm1-1,32);t3++) {
      for (t4=0;t4<=floord(32*t3+Nm1+30,32);t4++) {
        if ((t3 == 0) && (t4 == 0)) {
          for (t5=1;t5<=30;t5++) {
            Gim1 = G[t5 - 1];;
            Gip1 = G[t5 + 1];;
            Gi = G[t5];;
            for (t6=t5+1;t6<=31;t6++) {
              Gi[(-t5+t6)] = omega_over_four * (Gim1[(-t5+t6)] + Gip1[(-t5+t6)] + Gi[(-t5+t6) - 1] + Gi[(-t5+t6) + 1]) + one_minus_omega * Gi[(-t5+t6)];;
            }
          }
        }
        if ((t3 == 0) && (t4 == 0)) {
          Gim1 = G[31 - 1];;
          Gip1 = G[31 + 1];;
          Gi = G[31];;
        }
        if ((t3 == 1) && (t4 == 1)) {
          for (t5=32;t5<=Mm1-1;t5++) {
            Gi = G[t5];;
            for (t6=t5+1;t6<=63;t6++) {
              Gi[(-t5+t6)] = omega_over_four * (Gim1[(-t5+t6)] + Gip1[(-t5+t6)] + Gi[(-t5+t6) - 1] + Gi[(-t5+t6) + 1]) + one_minus_omega * Gi[(-t5+t6)];;
            }
          }
        }
        if ((t3 == 1) && (t4 == 0)) {
          for (t5=32;t5<=Mm1-1;t5++) {
            Gim1 = G[t5 - 1];;
            Gip1 = G[t5 + 1];;
          }
        }
        if (t3 <= t4-1) {
          for (t5=max(max(1,32*t3),32*t4-Nm1+1);t5<=min(Mm1-1,32*t3+31);t5++) {
            for (t6=32*t4;t6<=min(32*t4+31,t5+Nm1-1);t6++) {
              Gi[(-t5+t6)] = omega_over_four * (Gim1[(-t5+t6)] + Gip1[(-t5+t6)] + Gi[(-t5+t6) - 1] + Gi[(-t5+t6) + 1]) + one_minus_omega * Gi[(-t5+t6)];;
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
    ARRAY_PREPARATION_2D(array_0, M+1, L+1);
    ARRAY_PREPARATION_1D(array_1, L+1);
    ARRAY_PREPARATION_1D(array_2, L+1);
	ARRAY_PREPARATION_1D(array_3, L+1);

    SCImark_sor(array_0, array_1, array_2, array_3, N, M, L);

    print_array_1d(array_1, L+1);

    free_array_2d(array_0, M+1, L+1);
    free_array_1d(array_1, L+1);
    free_array_1d(array_2, L+1);
	free_array_1d(array_3, L+1);

    return 0;
}
