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
void NPB_bt5(double **ue, double ****forcing, int n_, int m_, int q_, int p_) {
    int i, j, k, m;
    double dssp = 0.1; // Assuming dssp is a scalar value used in the calculations  

    double time_start = omp_get_wtime();
#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++){
    for (int ii = 0; ii < n_; ii += 32) {
        for (int kk = 0; kk < q_; kk += 32) {
            for (int jj = 1 * 3; jj <= m_ - 3 * 1 - 1; jj += 32) {
                for (int mm = 0; mm < p_; mm += 32) {
                    for (i = ii; i < min(ii + 32, n_); i++) {
                        for (k = kk; k < min(kk + 32, q_); k++) {
                            for (j = jj; j < min(jj + 32, m_ - 3 * 1 - 1); j++) {
                                for (m = mm; m < min(mm + 32, p_); m++) {
                                    forcing[i][1][k][m] = forcing[i][1][k][m] - dssp *
                                        (5.0 * ue[1][m] - 4.0 * ue[1 + 1][m] + ue[1 + 2][m]);
                                    forcing[i][2][k][m] = forcing[i][2][k][m] - dssp *
                                        (-4.0 * ue[2 - 1][m] + 6.0 * ue[2][m] -
                                            4.0 * ue[2 + 1][m] + ue[2 + 2][m]);

                                    forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
                                        (ue[j - 2][m] - 4.0 * ue[j - 1][m] +
                                            6.0 * ue[j][m] - 4.0 * ue[j + 1][m] + ue[j + 2][m]);

                                    forcing[i][m_ - 3][k][m] = forcing[i][m_ - 3][k][m] - dssp *
                                        (ue[m_ - 3 - 2][m] - 4.0 * ue[m_ - 3 - 1][m] +
                                            6.0 * ue[m_ - 3][m] - 4.0 * ue[m_ - 3 + 1][m]);
                                    forcing[i][m_ - 2][k][m] = forcing[i][m_ - 2][k][m] - dssp *
                                        (ue[m_ - 2 - 2][m] - 4.0 * ue[m_ - 2 - 1][m] + 5.0 * ue[m_ - 2][m]);
                                }
                            }
                        }
                    }
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

int main() {

    ARRAY_PREPARATION_2D(array_0, M+1, P+1);
    ARRAY_PREPARATION_4D(array_1, N+1, M+1, Q+1, P+1);

    NPB_bt5(array_0, array_1, N, M, Q, P);
 
    print_array_4d(array_1, N+1, M+1, Q+1, P+1);

    free_array_2d(array_0, M+1, P+1);
    free_array_4d(array_1, N+1, M+1, Q+1, P+1);

    return 0;
}
