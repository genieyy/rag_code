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
#define N 80
#define M 100
#define L 5
/* end parameters define */

/* start kernel func */
void NPB_lu3(double ****d, double **tmat, int iend, int jend, int n)
{
	int i, j, m;
	int ist = 1;
	int jst = 2;

    double time_start = omp_get_wtime();
  int t1, t2, t3, t4, t5, t6;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
if ((ist <= iend) && (jst <= jend)) {
  for (t2=0;t2<=floord(ITERATIONS-1,32);t2++) {
    for (t3=0;t3<=n-1;t3++) {
      for (t4=32*t2;t4<=min(ITERATIONS-1,32*t2+31);t4++) {
        for (t5=ist;t5<=iend;t5++) {
          for (t6=jst;t6<=jend;t6++) {
            tmat[t3][0] = d[t5][t6][t3][0];;
            tmat[t3][1] = d[t5][t6][t3][1];;
            tmat[t3][2] = d[t5][t6][t3][2];;
            tmat[t3][3] = d[t5][t6][t3][3];;
            tmat[t3][4] = d[t5][t6][t3][4];;
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
	ARRAY_PREPARATION_4D(array_0, N + 1, M + 1, L + 1, 5);
	ARRAY_PREPARATION_2D(array_1, L + 1, 5);

	NPB_lu3(array_0, array_1, N, M, L);

	print_array_2d(array_1, L + 1, 5);

	free_array_4d(array_0, N + 1, M + 1, L + 1, 5);
	free_array_2d(array_1, L + 1, 5);

	return 0;
}
