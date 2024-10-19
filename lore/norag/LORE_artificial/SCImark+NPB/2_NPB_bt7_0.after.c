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
#define Q 40
#define P 5
/* end parameters define */

/* start kernel func */
void NPB_bt7(double **ue, double ****forcing, int n_, int m_, int q_, int p_) {
    int i, j, k, m;
    double dssp = 0.1; // Example value, replace with actual value  

    double time_start = omp_get_wtime();
#pragma scop
/**/

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (i = 0; i < n_; i++) {
        for (j = 0; j < m_; j++) {
            double dssp_times_5 = dssp * 5.0;
            double dssp_times_4 = dssp * 4.0;
            double dssp_times_6 = dssp * 6.0;

            for (m = 0; m < p_; m++) {
                forcing[i][j][1][m] -= dssp_times_5 * ue[1][m] - dssp_times_4 * (ue[2][m] + ue[3][m]);
                forcing[i][j][2][m] -= dssp_times_6 * ue[2][m] - dssp_times_4 * (ue[1][m] + ue[3][m]) + dssp * ue[4][m];
            }

            for (m = 0; m < p_; m++) {
                for (k = 3; k <= q_ - 4; k++) {
                    forcing[i][j][k][m] -= dssp * (ue[k - 2][m] - dssp_times_4 * (ue[k - 1][m] + ue[k + 1][m]) + dssp_times_6 * ue[k][m] + ue[k + 2][m]);
                }
            }

            for (m = 0; m < p_; m++) {
                forcing[i][j][q_ - 3][m] -= dssp * (ue[q_ - 5][m] - dssp_times_4 * (ue[q_ - 4][m] + ue[q_ - 2][m]) + dssp_times_6 * ue[q_ - 3][m]);
                forcing[i][j][q_ - 2][m] -= dssp * (ue[q_ - 4][m] - dssp_times_4 * ue[q_ - 3][m] + dssp_times_5 * ue[q_ - 2][m]);
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

    ARRAY_PREPARATION_2D(array_0, Q+1, P+1);
    ARRAY_PREPARATION_4D(array_1, N+1, M+1, Q+1, P+1);

    NPB_bt7(array_0, array_1, N, M, Q, P);
 
    print_array_4d(array_1, N+1, M+1, Q+1, P+1);

    free_array_2d(array_0, Q+1, P+1);
    free_array_4d(array_1, N+1, M+1, Q+1, P+1);

    return 0;
}