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

int t1, t2, t3, t4;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

lbp = 0;
ubp = floord(n_, 32);
#pragma omp parallel for private(lbv, ubv, t3, t4)
for (t1 = lbp; t1 <= ubp; t1++) {
    for (t2 = 0; t2 <= floord(32 * t1 + m_ - 1, 32); t2++) {
        for (t3 = max(32 * t1, 32 * t2); t3 <= min(n_ - 1, 32 * t1 + 31); t3++) {
            for (t4 = 32 * t2; t4 <= min(m_ - 1, 32 * t2 + 31); t4++) {
                for (m = 0; m < p_; m++) {
                    forcing[t3][t4][1][m] = forcing[t3][t4][1][m] - dssp *
                        (5.0 * ue[1][m] - 4.0 * ue[1 + 1][m] + ue[1 + 2][m]);
                    
                    forcing[t3][t4][2][m] = forcing[t3][t4][2][m] - dssp *
                        (-4.0 * ue[2 - 1][m] + 6.0 * ue[2][m] -
                            4.0 * ue[2 + 1][m] + ue[2 + 2][m]);
                }

                for (m = 0; m < p_; m++) {
                    for (k = 1 * 3; k <= q_ - 4; k++) {
                        forcing[t3][t4][k][m] = forcing[t3][t4][k][m] - dssp *
                            (ue[k - 2][m] - 4.0 * ue[k - 1][m] +
                                6.0 * ue[k][m] - 4.0 * ue[k + 1][m] + ue[k + 2][m]);
                    }
                }

                for (m = 0; m < p_; m++) {
                    forcing[t3][t4][q_ - 3][m] = forcing[t3][t4][q_ - 3][m] - dssp *
                        (ue[q_ - 3 - 2][m] - 4.0 * ue[q_ - 3 - 1][m] +
                            6.0 * ue[q_ - 3][m] - 4.0 * ue[q_ - 3 + 1][m]);
                    
                    forcing[t3][t4][q_ - 2][m] = forcing[t3][t4][q_ - 2][m] - dssp *
                        (ue[q_ - 2 - 2][m] - 4.0 * ue[q_ - 2 - 1][m] + 5.0 * ue[q_ - 2][m]);
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

    ARRAY_PREPARATION_2D(array_0, Q+1, P+1);
    ARRAY_PREPARATION_4D(array_1, N+1, M+1, Q+1, P+1);

    NPB_bt7(array_0, array_1, N, M, Q, P);
 
    print_array_4d(array_1, N+1, M+1, Q+1, P+1);

    free_array_2d(array_0, Q+1, P+1);
    free_array_4d(array_1, N+1, M+1, Q+1, P+1);

    return 0;
}