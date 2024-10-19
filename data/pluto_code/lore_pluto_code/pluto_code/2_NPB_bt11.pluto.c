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
void NPB_bt11(double ****u, double ***Pface, int nz, int mz, int q, int p)
{
    double xi, eta, zeta, Pxi, Peta, Pzeta;
    double dnxm1 = 0.1;  // Example values for dnxm1, dnym1, dnzm1  
    double dnym1 = 0.2;
    double dnzm1 = 0.3;

    double time_start = omp_get_wtime();
  int t1, t2, t3, t4, t5, t6;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
for (t1=0;t1<=ITERATIONS-1;t1++) {
  for (t2=0;t2<=nz-1;t2++) {
    zeta = (double)0 * dnzm1;;
    Pzeta = zeta * Pface[1][2][0] + (1.0 - zeta) * Pface[0][2][0];;
    eta = (double)0 * dnym1;;
    Peta = eta * Pface[1][1][0] + (1.0 - eta) * Pface[0][1][0];;
    xi = (double)t2 * dnxm1;;
    Pxi = xi * Pface[1][0][0] + (1.0 - xi) * Pface[0][0][0];;
    u[t2][0][0][0] = Pxi + Peta + Pzeta - Pxi * Peta - Pxi * Pzeta - Peta * Pzeta + Pxi * Peta * Pzeta;;
    for (t5=1;t5<=p-1;t5++) {
      Pzeta = zeta * Pface[1][2][t5] + (1.0 - zeta) * Pface[0][2][t5];;
      Peta = eta * Pface[1][1][t5] + (1.0 - eta) * Pface[0][1][t5];;
      Pxi = xi * Pface[1][0][t5] + (1.0 - xi) * Pface[0][0][t5];;
      u[t2][0][0][t5] = Pxi + Peta + Pzeta - Pxi * Peta - Pxi * Pzeta - Peta * Pzeta + Pxi * Peta * Pzeta;;
    }
    for (t4=1;t4<=q-1;t4++) {
      zeta = (double)t4 * dnzm1;;
      for (t5=0;t5<=p-1;t5++) {
        Pzeta = zeta * Pface[1][2][t5] + (1.0 - zeta) * Pface[0][2][t5];;
        Peta = eta * Pface[1][1][t5] + (1.0 - eta) * Pface[0][1][t5];;
        Pxi = xi * Pface[1][0][t5] + (1.0 - xi) * Pface[0][0][t5];;
        u[t2][0][t4][t5] = Pxi + Peta + Pzeta - Pxi * Peta - Pxi * Pzeta - Peta * Pzeta + Pxi * Peta * Pzeta;;
      }
    }
    for (t3=1;t3<=mz-1;t3++) {
      zeta = (double)0 * dnzm1;;
      Pzeta = zeta * Pface[1][2][0] + (1.0 - zeta) * Pface[0][2][0];;
      eta = (double)t3 * dnym1;;
      Peta = eta * Pface[1][1][0] + (1.0 - eta) * Pface[0][1][0];;
      Pxi = xi * Pface[1][0][0] + (1.0 - xi) * Pface[0][0][0];;
      u[t2][t3][0][0] = Pxi + Peta + Pzeta - Pxi * Peta - Pxi * Pzeta - Peta * Pzeta + Pxi * Peta * Pzeta;;
      for (t5=1;t5<=p-1;t5++) {
        Pzeta = zeta * Pface[1][2][t5] + (1.0 - zeta) * Pface[0][2][t5];;
        Peta = eta * Pface[1][1][t5] + (1.0 - eta) * Pface[0][1][t5];;
        Pxi = xi * Pface[1][0][t5] + (1.0 - xi) * Pface[0][0][t5];;
        u[t2][t3][0][t5] = Pxi + Peta + Pzeta - Pxi * Peta - Pxi * Pzeta - Peta * Pzeta + Pxi * Peta * Pzeta;;
      }
      for (t4=1;t4<=q-1;t4++) {
        zeta = (double)t4 * dnzm1;;
        for (t5=0;t5<=p-1;t5++) {
          Pzeta = zeta * Pface[1][2][t5] + (1.0 - zeta) * Pface[0][2][t5];;
          Peta = eta * Pface[1][1][t5] + (1.0 - eta) * Pface[0][1][t5];;
          Pxi = xi * Pface[1][0][t5] + (1.0 - xi) * Pface[0][0][t5];;
          u[t2][t3][t4][t5] = Pxi + Peta + Pzeta - Pxi * Peta - Pxi * Pzeta - Peta * Pzeta + Pxi * Peta * Pzeta;;
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
    ARRAY_PREPARATION_4D(array_0, N + 1, M + 1, L + 1, P + 1);
    ARRAY_PREPARATION_3D(array_1, 2, 3, P + 1);

    NPB_bt11(array_0, array_1, N, M, L, P);

    print_array_4d(array_0, N + 1, M + 1, L + 1, P + 1);

    free_array_4d(array_0, N + 1, M + 1, L + 1, P + 1);
    free_array_3d(array_1, 2, 3, P + 1);

    return 0;
}
