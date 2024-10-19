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
void NPB_mg3(double ***r, double *r1, double * r2, double ***u, int n3, int n2, int n1)
{
	int i1, i2, i3;
	double c[3] = {0.2, 0.3, 0.4};

    double time_start = omp_get_wtime();
  int t1, t2, t3, t4, t5, t6, t7;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
for (t1=0;t1<=ITERATIONS-1;t1++) {
  for (t2=1;t2<=n3-2;t2++) {
    for (t3=0;t3<=floord(3*n2+n1-7,32);t3++) {
      lbp=max(ceild(2*t3,3),ceild(32*t3-n2+2,32));
      ubp=min(min(floord(2*n2+n1-5,32),floord(64*t3+n1+61,96)),t3);
#pragma omp parallel for private(lbv,ubv,t5,t6,t7)
      for (t4=lbp;t4<=ubp;t4++) {
        for (t5=max(max(1,ceild(32*t4-n1+1,2)),32*t3-32*t4);t5<=min(min(n2-2,16*t4+14),32*t3-32*t4+31);t5++) {
          for (t6=max(32*t4,2*t5);t6<=2*t5+1;t6++) {
            r2[(-2*t5+t6)] = r[t2 - 1][t5 - 1][(-2*t5+t6)] + r[t2 - 1][t5 + 1][(-2*t5+t6)] + r[t2 + 1][t5 - 1][(-2*t5+t6)] + r[t2 + 1][t5 + 1][(-2*t5+t6)];;
            r1[(-2*t5+t6)] = r[t2][t5 - 1][(-2*t5+t6)] + r[t2][t5 + 1][(-2*t5+t6)] + r[t2 - 1][t5][(-2*t5+t6)] + r[t2 + 1][t5][(-2*t5+t6)];;
          }
          for (t6=max(32*t4,2*t5+2);t6<=min(32*t4+31,2*t5+n1-1);t6++) {
            r2[(-2*t5+t6)] = r[t2 - 1][t5 - 1][(-2*t5+t6)] + r[t2 - 1][t5 + 1][(-2*t5+t6)] + r[t2 + 1][t5 - 1][(-2*t5+t6)] + r[t2 + 1][t5 + 1][(-2*t5+t6)];;
            r1[(-2*t5+t6)] = r[t2][t5 - 1][(-2*t5+t6)] + r[t2][t5 + 1][(-2*t5+t6)] + r[t2 - 1][t5][(-2*t5+t6)] + r[t2 + 1][t5][(-2*t5+t6)];;
            u[t2][t5][(-2*t5+t6-1)] = u[t2][t5][(-2*t5+t6-1)] + c[0] * r[t2][t5][(-2*t5+t6-1)] + c[1] * (r[t2][t5][(-2*t5+t6-1) - 1] + r[t2][t5][(-2*t5+t6-1) + 1] + r1[(-2*t5+t6-1)]) + c[2] * (r2[(-2*t5+t6-1)] + r1[(-2*t5+t6-1) - 1] + r1[(-2*t5+t6-1) + 1]);;
          }
        }
        if ((t3 >= ceild(3*t4-1,2)) && (t4 <= floord(n2-17,16))) {
          for (t6=32*t4+30;t6<=32*t4+31;t6++) {
            r2[(-32*t4+t6-30)] = r[t2 - 1][(16*t4+15) - 1][(-32*t4+t6-30)] + r[t2 - 1][(16*t4+15) + 1][(-32*t4+t6-30)] + r[t2 + 1][(16*t4+15) - 1][(-32*t4+t6-30)] + r[t2 + 1][(16*t4+15) + 1][(-32*t4+t6-30)];;
            r1[(-32*t4+t6-30)] = r[t2][(16*t4+15) - 1][(-32*t4+t6-30)] + r[t2][(16*t4+15) + 1][(-32*t4+t6-30)] + r[t2 - 1][(16*t4+15)][(-32*t4+t6-30)] + r[t2 + 1][(16*t4+15)][(-32*t4+t6-30)];;
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
	ARRAY_PREPARATION_3D(array_0, N + 1, M + 1, L + 1);
	ARRAY_PREPARATION_1D(array_1, L+1);
	ARRAY_PREPARATION_1D(array_2, L+1);
	ARRAY_PREPARATION_3D(array_3, N + 1, M + 1, L + 1);

	NPB_mg3(array_0, array_1, array_2, array_3, N, M, L);

	print_array_3d(array_3, N + 1, M + 1, L + 1);

	free_array_3d(array_0, N + 1, M + 1, L + 1);
	free_array_1d(array_1, L+1);
	free_array_1d(array_2, L+1);
	free_array_3d(array_3, N + 1, M + 1, L + 1);

	return 0;
}
