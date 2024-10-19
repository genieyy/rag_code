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
#define N 20
#define M 30
#define L 40
#define P 5
/* end parameters define */

/* start kernel func */
void NPB_lu(double ****v, double ****ldz, int iend, int jend, int nz, int n)
{
  int i, j, k, m;
  int ist = 1;
  int jst = 2;
  double omega = 0.1;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation of the Optimization:
1. **Loop Unrolling and Temporary Variable Usage**:
   - The inner loop over `m` is unrolled to calculate the sum of products for each `m` index into a temporary array `temp`. This reduces the number of redundant calculations by storing the intermediate results in `temp`.
   - This approach minimizes the number of accesses to the `v` and `ldz` arrays, which can be costly in terms of memory access latency.

2. **Reduced Memory Accesses**:
   - By storing the intermediate results in `temp`, we avoid recalculating the same values multiple times within the loop over `m`. This reduces the overall number of memory accesses, which can significantly improve performance, especially if the arrays `v` and `ldz` are large.

3. **Improved Locality of Reference**:
   - The temporary array `temp` is used to store the results of the inner loop, which improves the locality of reference. This means that the CPU cache is more likely to hold the relevant data, reducing cache misses and improving performance.

This optimization leverages the principles of loop unrolling and temporary variable usage to enhance the performance of the original code.*/

for (int iter = 0; iter < ITERATIONS; iter++){
  for (i = ist; i <= iend; i++)
  {
    for (j = jst; j <= jend; j++)
    {
      for (k = 1; k < nz - 1; k++)
      {
        double temp[5];
        for (int idx = 0; idx < 5; idx++) {
          temp[idx] = ldz[i][j][0][idx] * v[i][j][k - 1][0] +
                      ldz[i][j][1][idx] * v[i][j][k - 1][1] +
                      ldz[i][j][2][idx] * v[i][j][k - 1][2] +
                      ldz[i][j][3][idx] * v[i][j][k - 1][3] +
                      ldz[i][j][4][idx] * v[i][j][k - 1][4];
        }
        for (m = 0; m < n; m++)
        {
          v[i][j][k][m] -= omega * temp[m];
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
  ARRAY_PREPARATION_4D(array_0, N + 1, M + 1, L + 1, P + 1);
  ARRAY_PREPARATION_4D(array_1, N + 1, M + 1, L + 1, 5);

  NPB_lu(array_0, array_1, N, M, L, P);

  print_array_4d(array_0, N + 1, M + 1, L + 1, P + 1);

  free_array_4d(array_0, N + 1, M + 1, L + 1, P + 1);
  free_array_4d(array_1, N + 1, M + 1, L + 1, 5);

  return 0;
}