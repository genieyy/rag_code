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

int t1, t2, t3, t4, t5, t6;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

for (int iter = 0; iter < ITERATIONS; iter++) {
    lbp = 0;
    ubp = floord(n_, 32);
#pragma omp parallel for private(lbv, ubv, t3, t4, t5, t6)
    for (t1 = lbp; t1 <= ubp; t1++) {
        for (t2 = 0; t2 <= floord(m_, 32); t2++) {
            for (t3 = 0; t3 <= floord(p_, 32); t3++) {
                for (t4 = max(32 * t1, 32 * t2); t4 <= min(n_, 32 * t1 + 31); t4++) {
                    for (t5 = 32 * t2; t5 <= min(m_, 32 * t2 + 31); t5++) {
                        for (t6 = 32 * t3; t6 <= min(p_, 32 * t3 + 31); t6++) {
                            forcing[t4][t5][1][t6] = forcing[t4][t5][1][t6] - dssp *
                                (5.0 * ue[1][t6] - 4.0 * ue[1 + 1][t6] + ue[1 + 2][t6]);
                            
                            forcing[t4][t5][2][t6] = forcing[t4][t5][2][t6] - dssp *
                                (-4.0 * ue[2 - 1][t6] + 6.0 * ue[2][t6] -
                                    4.0 * ue[2 + 1][t6] + ue[2 + 2][t6]);
                        }
                    }
                }
                for (t4 = max(32 * t1, 32 * t2); t4 <= min(n_, 32 * t1 + 31); t4++) {
                    for (t5 = 32 * t2; t5 <= min(m_, 32 * t2 + 31); t5++) {
                        for (t6 = 32 * t3; t6 <= min(p_, 32 * t3 + 31); t6++) {
                            for (int k = 1 * 3; k <= q_ - 4; k++) {
                                forcing[t4][t5][k][t6] = forcing[t4][t5][k][t6] - dssp *
                                    (ue[k - 2][t6] - 4.0 * ue[k - 1][t6] +
                                        6.0 * ue[k][t6] - 4.0 * ue[k + 1][t6] + ue[k + 2][t6]);
                            }
                        }
                    }
                }
                for (t4 = max(32 * t1, 32 * t2); t4 <= min(n_, 32 * t1 + 31); t4++) {
                    for (t5 = 32 * t2; t5 <= min(m_, 32 * t2 + 31); t5++) {
                        for (t6 = 32 * t3; t6 <= min(p_, 32 * t3 + 31); t6++) {
                            forcing[t4][t5][q_ - 3][t6] = forcing[t4][t5][q_ - 3][t6] - dssp *
                                (ue[q_ - 3 - 2][t6] - 4.0 * ue[q_ - 3 - 1][t6] +
                                    6.0 * ue[q_ - 3][t6] - 4.0 * ue[q_ - 3 + 1][t6]);
                            
                            forcing[t4][t5][q_ - 2][t6] = forcing[t4][t5][q_ - 2][t6] - dssp *
                                (ue[q_ - 2 - 2][t6] - 4.0 * ue[q_ - 2 - 1][t6] + 5.0 * ue[q_ - 2][t6]);
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

    ARRAY_PREPARATION_2D(array_0, Q+1, P+1);
    ARRAY_PREPARATION_4D(array_1, N+1, M+1, Q+1, P+1);

    NPB_bt7(array_0, array_1, N, M, Q, P);
 
    print_array_4d(array_1, N+1, M+1, Q+1, P+1);

    free_array_2d(array_0, Q+1, P+1);
    free_array_4d(array_1, N+1, M+1, Q+1, P+1);

    return 0;
}