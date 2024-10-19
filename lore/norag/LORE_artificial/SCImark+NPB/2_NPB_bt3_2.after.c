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
#define N 10
#define M 20
#define Q 30
#define P 5
/* end parameters define */

/* start kernel func */
void NPB_bt3(double *q, double *cuf, double **ue, double **buf, double ****forcing, double *dtemp, int n_, int m_, int q_, int p_)
{
    double dnym1, dnzm1, dnxm1;
    double eta, zeta, xi, dtpp;
    int j, k, i, m;
    double tx2, dx1tx1, xxcon1, dx2tx1, xxcon2, dx3tx1, xxcon3, xxcon4, xxcon5, dx5tx1;
    double c1, c2;


    tx2 = 0.1; dx1tx1 = 0.2; xxcon1 = 0.3; dx2tx1 = 0.4; xxcon2 = 0.5;
    dx3tx1 = 0.6; xxcon3 = 0.7; xxcon4 = 0.8; xxcon5 = 0.9; dx5tx1 = 1.0;
    double dx4tx1 = 0.8;
    c1 = 1.1; c2 = 1.2;

    double time_start = omp_get_wtime();
#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++) {
    for (j = 1; j < n_ - 1; j++) {
        eta = (double)j * dnym1;
        for (k = 1; k < m_ - 1; k++) {
            zeta = (double)k * dnzm1;
            for (i = 0; i < q_; i++) {
                xi = (double)i * dnxm1;
                for (m = 0; m < p_; m++) {
                    ue[i][m] = dtemp[m];
                }
                dtpp = 1.0 / dtemp[0];
                for (m = 1; m <= p_ - 1; m++) {
                    buf[i][m] = dtpp * dtemp[m];
                }
                cuf[i] = buf[i][1] * buf[i][1];
                buf[i][0] = cuf[i] + buf[i][2] * buf[i][2] + buf[i][3] * buf[i][3];
                q[i] = 0.5 * (buf[i][1] * ue[i][1] + buf[i][2] * ue[i][2] + buf[i][3] * ue[i][3]);
            }
            for (i = 1; i < q_ - 1; i += 2) {
                forcing[i][j][k][0] = forcing[i][j][k][0] - tx2 * (ue[i+1][1] - ue[i-1][1]) + dx1tx1 * (ue[i+1][0] - 2.0 * ue[i][0] + ue[i-1][0]);
                forcing[i][j][k][1] = forcing[i][j][k][1] - tx2 * ((ue[i+1][1] * buf[i+1][1] + c2 * (ue[i+1][4] - q[i+1])) - (ue[i-1][1] * buf[i-1][1] + c2 * (ue[i-1][4] - q[i-1]))) + xxcon1 * (buf[i+1][1] - 2.0 * buf[i][1] + buf[i-1][1]) + dx2tx1 * (ue[i+1][1] - 2.0 * ue[i][1] + ue[i-1][1]);
                forcing[i][j][k][2] = forcing[i][j][k][2] - tx2 * (ue[i+1][2] * buf[i+1][1] - ue[i-1][2] * buf[i-1][1]) + xxcon2 * (buf[i+1][2] - 2.0 * buf[i][2] + buf[i-1][2]) + dx3tx1 * (ue[i+1][2] - 2.0 * ue[i][2] + ue[i-1][2]);
                forcing[i][j][k][3] = forcing[i][j][k][3] - tx2 * (ue[i+1][3] * buf[i+1][1] - ue[i-1][3] * buf[i-1][1]) + xxcon2 * (buf[i+1][3] - 2.0 * buf[i][3] + buf[i-1][3]) + dx4tx1 * (ue[i+1][3] - 2.0 * ue[i][3] + ue[i-1][3]);
                forcing[i][j][k][4] = forcing[i][j][k][4] - tx2 * (buf[i+1][1] * (c1 * ue[i+1][4] - c2 * q[i+1]) - buf[i-1][1] * (c1 * ue[i-1][4] - c2 * q[i-1])) + 0.5 * xxcon3 * (buf[i+1][0] - 2.0 * buf[i][0] + buf[i-1][0]) + xxcon4 * (cuf[i+1] - 2.0 * cuf[i] + cuf[i-1]) + xxcon5 * (buf[i+1][4] - 2.0 * buf[i][4] + buf[i-1][4]) + dx5tx1 * (ue[i+1][4] - 2.0 * ue[i][4] + ue[i-1][4]);

                forcing[i+1][j][k][0] = forcing[i+1][j][k][0] - tx2 * (ue[i+2][1] - ue[i][1]) + dx1tx1 * (ue[i+2][0] - 2.0 * ue[i+1][0] + ue[i][0]);
                forcing[i+1][j][k][1] = forcing[i+1][j][k][1] - tx2 * ((ue[i+2][1] * buf[i+2][1] + c2 * (ue[i+2][4] - q[i+2])) - (ue[i][1] * buf[i][1] + c2 * (ue[i][4] - q[i]))) + xxcon1 * (buf[i+2][1] - 2.0 * buf[i+1][1] + buf[i][1]) + dx2tx1 * (ue[i+2][1] - 2.0 * ue[i+1][1] + ue[i][1]);
                forcing[i+1][j][k][2] = forcing[i+1][j][k][2] - tx2 * (ue[i+2][2] * buf[i+2][1] - ue[i][2] * buf[i][1]) + xxcon2 * (buf[i+2][2] - 2.0 * buf[i+1][2] + buf[i][2]) + dx3tx1 * (ue[i+2][2] - 2.0 * ue[i+1][2] + ue[i][2]);
                forcing[i+1][j][k][3] = forcing[i+1][j][k][3] - tx2 * (ue[i+2][3] * buf[i+2][1] - ue[i][3] * buf[i][1]) + xxcon2 * (buf[i+2][3] - 2.0 * buf[i+1][3] + buf[i][3]) + dx4tx1 * (ue[i+2][3] - 2.0 * ue[i+1][3] + ue[i][3]);
                forcing[i+1][j][k][4] = forcing[i+1][j][k][4] - tx2 * (buf[i+2][1] * (c1 * ue[i+2][4] - c2 * q[i+2]) - buf[i][1] * (c1 * ue[i][4] - c2 * q[i]))) + 0.5 * xxcon3 * (buf[i+2][0] - 2.0 * buf[i+1][0] + buf[i][0]) + xxcon4 * (cuf[i+2] - 2.0 * cuf[i+1] + cuf[i]) + xxcon5 * (buf[i+2][4] - 2.0 * buf[i+1][4] + buf[i][4]) + dx5tx1 * (ue[i+2][4] - 2.0 * ue[i+1][4] + ue[i][4]);
            }
        }
    }
}
#pragma endscop
    double time_end = omp_get_wtime();
    printf("%f\n", time_end - time_start);
}
/* end kernel func */

int main() {

    ARRAY_PREPARATION_1D(array_0, Q+1);
    ARRAY_PREPARATION_1D(array_1, Q+1);
    ARRAY_PREPARATION_2D(array_2, Q+1, P+1);
    ARRAY_PREPARATION_2D(array_3, Q+1, P+1);
    ARRAY_PREPARATION_4D(array_4, Q+1, N+1, M+1, P+1);
    ARRAY_PREPARATION_1D(array_5, P+1);

    NPB_bt3(array_0, array_1, array_2, array_3, array_4, array_5, N, M, Q, P);

    print_array_4d(array_4, Q+1, N+1, M+1, P+1);

    free_array_1d(array_0, Q+1);
    free_array_1d(array_1, Q+1);
    free_array_2d(array_2, Q+1, P+1);
    free_array_2d(array_3, Q+1, P+1);
    free_array_4d(array_4, Q+1, N+1, M+1, P+1);
    free_array_1d(array_5, P+1);

    return 0;
}