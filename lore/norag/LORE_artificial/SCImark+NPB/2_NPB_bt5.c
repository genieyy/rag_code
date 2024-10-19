#include <stdio.h>  
#include <stdlib.h>  
#include <string.h>
#include <math.h>
#include "lore.h"

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

#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++){
    for (i = 0; i < n_; i++) {
        for (k = 0; k < q_; k++) {
            for (m = 0; m < p_; m++) {
                forcing[i][1][k][m] = forcing[i][1][k][m] - dssp *
                    (5.0 * ue[1][m] - 4.0 * ue[1 + 1][m] + ue[1 + 2][m]);
            
                forcing[i][2][k][m] = forcing[i][2][k][m] - dssp *
                    (-4.0 * ue[2 - 1][m] + 6.0 * ue[2][m] -
                        4.0 * ue[2 + 1][m] + ue[2 + 2][m]);
            }

            for (m = 0; m < p_; m++) {
                for (j = 1 * 3; j <= m_ - 3 * 1 - 1; j++) {
                    forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
                        (ue[j - 2][m] - 4.0 * ue[j - 1][m] +
                            6.0 * ue[j][m] - 4.0 * ue[j + 1][m] + ue[j + 2][m]);
                }
            }
            
            for (m = 0; m < p_; m++) {
                forcing[i][m_ - 3][k][m] = forcing[i][m_ - 3][k][m] - dssp *
                    (ue[m_ - 3 - 2][m] - 4.0 * ue[m_ - 3 - 1][m] +
                        6.0 * ue[m_ - 3][m] - 4.0 * ue[m_ - 3 + 1][m]);
            
                forcing[i][m_ - 2][k][m] = forcing[i][m_ - 2][k][m] - dssp *
                    (ue[m_ - 2 - 2][m] - 4.0 * ue[m_ - 2 - 1][m] + 5.0 * ue[m_ - 2][m]);
            }
        }
    }
}
#pragma endscop
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
