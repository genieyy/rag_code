#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lore.h"

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

#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++){
  for (i = ist; i <= iend; i++)
  {
    for (j = jst; j <= jend; j++)
    {
      for (k = 1; k < nz - 1; k++)
      {
        for (m = 0; m < n; m++)
        {
          v[i][j][k][m] = v[i][j][k][m] - omega * (ldz[i][j][m][0] * v[i][j][k - 1][0] + ldz[i][j][m][1] * v[i][j][k - 1][1] + ldz[i][j][m][2] * v[i][j][k - 1][2] + ldz[i][j][m][3] * v[i][j][k - 1][3] + ldz[i][j][m][4] * v[i][j][k - 1][4]);
        }
      }
    }
  }
}
#pragma endscop
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