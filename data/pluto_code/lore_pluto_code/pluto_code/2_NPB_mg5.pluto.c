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
void NPB_mg5(double ***u, double ***z, int mm3, int mm2, int mm1)
{
	int i1, i2, i3;
	int d1 = 2;
	int d2 = 3;
	int c1 = 1;
	int c2 = 1;
	int c3 = 1;

    double time_start = omp_get_wtime();
  int t1, t2, t3, t4, t5, t6, t7, t8, t9;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
lbp=0;
ubp=floord(mm3-1,32);
#pragma omp parallel for private(lbv,ubv,t2,t3,t4,t5,t6,t7,t8,t9)
for (t1=lbp;t1<=ubp;t1++) {
  for (t2=0;t2<=floord(ITERATIONS-1,32);t2++) {
    if ((d1 <= mm1-1) && (d2 <= mm2-1)) {
      for (t3=max(1,32*t1);t3<=min(mm3-1,32*t1+31);t3++) {
        for (t4=32*t2;t4<=min(ITERATIONS-1,32*t2+31);t4++) {
          for (t6=d2;t6<=mm2-1;t6++) {
            lbv=d1;
            ubv=mm1-1;
#pragma ivdep
#pragma vector always
            for (t9=lbv;t9<=ubv;t9++) {
              u[2 * t3 - c3 - 1][2 * t6 - d2 - 1][2 * t9 - d1 - 1] = u[2 * t3 - c3 - 1][2 * t6 - d2 - 1][2 * t9 - d1 - 1] + 0.5 * (z[t3][t6 - 1][t9 - 1] + z[t3 - 1][t6 - 1][t9 - 1]);;
            }
            lbv=1;
            ubv=mm1-1;
#pragma ivdep
#pragma vector always
            for (t9=lbv;t9<=ubv;t9++) {
              u[2 * t3 - c3 - 1][2 * t6 - d2 - 1][2 * t9 - c1 - 1] = u[2 * t3 - c3 - 1][2 * t6 - d2 - 1][2 * t9 - c1 - 1] + 0.25 * (z[t3][t6 - 1][t9] + z[t3][t6 - 1][t9 - 1] + z[t3 - 1][t6 - 1][t9] + z[t3 - 1][t6 - 1][t9 - 1]);;
            }
          }
          for (t6=1;t6<=mm2-1;t6++) {
            lbv=d1;
            ubv=mm1-1;
#pragma ivdep
#pragma vector always
            for (t9=lbv;t9<=ubv;t9++) {
              u[2 * t3 - c3 - 1][2 * t6 - c2 - 1][2 * t9 - d1 - 1] = u[2 * t3 - c3 - 1][2 * t6 - c2 - 1][2 * t9 - d1 - 1] + 0.25 * (z[t3][t6][t9 - 1] + z[t3][t6 - 1][t9 - 1] + z[t3 - 1][t6][t9 - 1] + z[t3 - 1][t6 - 1][t9 - 1]);;
            }
            lbv=1;
            ubv=mm1-1;
#pragma ivdep
#pragma vector always
            for (t9=lbv;t9<=ubv;t9++) {
              u[2 * t3 - c3 - 1][2 * t6 - c2 - 1][2 * t9 - c1 - 1] = u[2 * t3 - c3 - 1][2 * t6 - c2 - 1][2 * t9 - c1 - 1] + 0.125 * (z[t3][t6][t9] + z[t3][t6 - 1][t9] + z[t3][t6][t9 - 1] + z[t3][t6 - 1][t9 - 1] + z[t3 - 1][t6][t9] + z[t3 - 1][t6 - 1][t9] + z[t3 - 1][t6][t9 - 1] + z[t3 - 1][t6 - 1][t9 - 1]);;
            }
          }
        }
      }
    }
    if ((d1 >= mm1) && (d2 <= mm2-1)) {
      for (t3=max(1,32*t1);t3<=min(mm3-1,32*t1+31);t3++) {
        for (t4=32*t2;t4<=min(ITERATIONS-1,32*t2+31);t4++) {
          for (t6=d2;t6<=mm2-1;t6++) {
            lbv=1;
            ubv=mm1-1;
#pragma ivdep
#pragma vector always
            for (t9=lbv;t9<=ubv;t9++) {
              u[2 * t3 - c3 - 1][2 * t6 - d2 - 1][2 * t9 - c1 - 1] = u[2 * t3 - c3 - 1][2 * t6 - d2 - 1][2 * t9 - c1 - 1] + 0.25 * (z[t3][t6 - 1][t9] + z[t3][t6 - 1][t9 - 1] + z[t3 - 1][t6 - 1][t9] + z[t3 - 1][t6 - 1][t9 - 1]);;
            }
          }
          for (t6=1;t6<=mm2-1;t6++) {
            lbv=1;
            ubv=mm1-1;
#pragma ivdep
#pragma vector always
            for (t9=lbv;t9<=ubv;t9++) {
              u[2 * t3 - c3 - 1][2 * t6 - c2 - 1][2 * t9 - c1 - 1] = u[2 * t3 - c3 - 1][2 * t6 - c2 - 1][2 * t9 - c1 - 1] + 0.125 * (z[t3][t6][t9] + z[t3][t6 - 1][t9] + z[t3][t6][t9 - 1] + z[t3][t6 - 1][t9 - 1] + z[t3 - 1][t6][t9] + z[t3 - 1][t6 - 1][t9] + z[t3 - 1][t6][t9 - 1] + z[t3 - 1][t6 - 1][t9 - 1]);;
            }
          }
        }
      }
    }
    if ((d1 <= mm1-1) && (d2 >= mm2)) {
      for (t3=max(1,32*t1);t3<=min(mm3-1,32*t1+31);t3++) {
        for (t4=32*t2;t4<=min(ITERATIONS-1,32*t2+31);t4++) {
          for (t6=1;t6<=mm2-1;t6++) {
            lbv=d1;
            ubv=mm1-1;
#pragma ivdep
#pragma vector always
            for (t9=lbv;t9<=ubv;t9++) {
              u[2 * t3 - c3 - 1][2 * t6 - c2 - 1][2 * t9 - d1 - 1] = u[2 * t3 - c3 - 1][2 * t6 - c2 - 1][2 * t9 - d1 - 1] + 0.25 * (z[t3][t6][t9 - 1] + z[t3][t6 - 1][t9 - 1] + z[t3 - 1][t6][t9 - 1] + z[t3 - 1][t6 - 1][t9 - 1]);;
            }
            lbv=1;
            ubv=mm1-1;
#pragma ivdep
#pragma vector always
            for (t9=lbv;t9<=ubv;t9++) {
              u[2 * t3 - c3 - 1][2 * t6 - c2 - 1][2 * t9 - c1 - 1] = u[2 * t3 - c3 - 1][2 * t6 - c2 - 1][2 * t9 - c1 - 1] + 0.125 * (z[t3][t6][t9] + z[t3][t6 - 1][t9] + z[t3][t6][t9 - 1] + z[t3][t6 - 1][t9 - 1] + z[t3 - 1][t6][t9] + z[t3 - 1][t6 - 1][t9] + z[t3 - 1][t6][t9 - 1] + z[t3 - 1][t6 - 1][t9 - 1]);;
            }
          }
        }
      }
    }
    if ((d1 >= mm1) && (d2 >= mm2)) {
      for (t3=max(1,32*t1);t3<=min(mm3-1,32*t1+31);t3++) {
        for (t4=32*t2;t4<=min(ITERATIONS-1,32*t2+31);t4++) {
          for (t6=1;t6<=mm2-1;t6++) {
            lbv=1;
            ubv=mm1-1;
#pragma ivdep
#pragma vector always
            for (t9=lbv;t9<=ubv;t9++) {
              u[2 * t3 - c3 - 1][2 * t6 - c2 - 1][2 * t9 - c1 - 1] = u[2 * t3 - c3 - 1][2 * t6 - c2 - 1][2 * t9 - c1 - 1] + 0.125 * (z[t3][t6][t9] + z[t3][t6 - 1][t9] + z[t3][t6][t9 - 1] + z[t3][t6 - 1][t9 - 1] + z[t3 - 1][t6][t9] + z[t3 - 1][t6 - 1][t9] + z[t3 - 1][t6][t9 - 1] + z[t3 - 1][t6 - 1][t9 - 1]);;
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
	ARRAY_PREPARATION_3D(array_0, 2*N+1, 2*M+1, 2*L+1);
	ARRAY_PREPARATION_3D(array_1, 2*N+1, 2*M+1, 2*L+1);

	NPB_mg5(array_0, array_1, N, M, L);

	print_array_3d(array_0, 2*N+1, 2*M+1, 2*L+1);

	free_array_3d(array_0, 2*N+1, 2*M+1, 2*L+1);
	free_array_3d(array_1, 2*N+1, 2*M+1, 2*L+1);

	return 0;
}
