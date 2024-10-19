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
void NPB_bt6(double **ue, double ****forcing, int n_, int m_, int q_, int p_) {
    int i, j, k, m;
    double dssp = 0.1; // Example value, replace with actual value  
    
    double time_start = omp_get_wtime();
#pragma scop
/**/

for (int iter = 0; iter < ITERATIONS; iter++){
    for (i = 0; i < n_; i++) {
        for (k = 0; k < q_; k++) {
            // Precompute common terms
            double dssp_5 = 5.0 * dssp;
            double dssp_4 = 4.0 * dssp;
            double dssp_6 = 6.0 * dssp;

            // Optimized inner loops
            for (m = 0; m < p_; m++) {
                double ue_1_m = ue[1][m];
                double ue_2_m = ue[2][m];
                double ue_3_m = ue[3][m];
                double ue_4_m = ue[4][m];

                forcing[i][1][k][m] -= dssp_5 * ue_1_m - dssp_4 * ue_2_m + dssp * ue_3_m;
                forcing[i][2][k][m] -= -dssp_4 * ue_1_m + dssp_6 * ue_2_m - dssp_4 * ue_3_m + dssp * ue_4_m;
            }

            for (m = 0; m < p_; m++) {
                for (j = 3; j <= m_ - 4; j++) {
                    double ue_jm2 = ue[j - 2][m];
                    double ue_jm1 = ue[j - 1][m];
                    double ue_j = ue[j][m];
                    double ue_jp1 = ue[j + 1][m];
                    double ue_jp2 = ue[j + 2][m];

                    forcing[i][j][k][m] -= dssp * (ue_jm2 - dssp_4 * ue_jm1 + dssp_6 * ue_j - dssp_4 * ue_jp1 + ue_jp2);
                }
            }

            for (m = 0; m < p_; m++) {
                double ue_m_3_m2 = ue[m_ - 5][m];
                double ue_m_3_m1 = ue[m_ - 4][m];
                double ue_m_3 = ue[m_ - 3][m];
                double ue_m_2 = ue[m_ - 2][m];
                double ue_m_1 = ue[m_ - 1][m];

                forcing[i][m_ - 3][k][m] -= dssp * (ue_m_3_m2 - dssp_4 * ue_m_3_m1 + dssp_6 * ue_m_3 - dssp_4 * ue_m_2);
                forcing[i][m_ - 2][k][m] -= dssp * (ue_m_3_m2 - dssp_4 * ue_m_3_m1 + dssp_5 * ue_m_2);
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

    ARRAY_PREPARATION_2D(array_0, M+1, P+1);
    ARRAY_PREPARATION_4D(array_1, N+1, M+1, Q+1, P+1);

    NPB_bt6(array_0, array_1, N, M, Q, P);
 
    print_array_4d(array_1, N+1, M+1, Q+1, P+1);

    free_array_2d(array_0, M+1, P+1);
    free_array_4d(array_1, N+1, M+1, Q+1, P+1);

    return 0;
}