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
/**/

int lb, ub, lbp, ubp;
register int lbv, ubv;

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int t1 = 0; t1 <= floord(n_ - 1, 32); t1++) {
        for (int t2 = 0; t2 <= floord(q_ - 1, 32); t2++) {
            for (int t3 = 0; t3 <= floord(p_ - 1, 32); t3++) {
                for (int t4 = max(0, 32 * t1); t4 <= min(n_ - 1, 32 * t1 + 31); t4++) {
                    for (int t5 = max(0, 32 * t2); t5 <= min(q_ - 1, 32 * t2 + 31); t5++) {
                        lbv = max(0, 32 * t3);
                        ubv = min(p_ - 1, 32 * t3 + 31);
#pragma ivdep
#pragma vector always
                        for (int t6 = lbv; t6 <= ubv; t6++) {
                            forcing[t4][1][t5][t6] = forcing[t4][1][t5][t6] - dssp *
                                (5.0 * ue[1][t6] - 4.0 * ue[1 + 1][t6] + ue[1 + 2][t6]);
                            
                            forcing[t4][2][t5][t6] = forcing[t4][2][t5][t6] - dssp *
                                (-4.0 * ue[2 - 1][t6] + 6.0 * ue[2][t6] -
                                    4.0 * ue[2 + 1][t6] + ue[2 + 2][t6]);
                        }

                        for (int t6 = lbv; t6 <= ubv; t6++) {
                            for (int t7 = 1 * 3; t7 <= m_ - 3 * 1 - 1; t7++) {
                                forcing[t4][t7][t5][t6] = forcing[t4][t7][t5][t6] - dssp *
                                    (ue[t7 - 2][t6] - 4.0 * ue[t7 - 1][t6] +
                                        6.0 * ue[t7][t6] - 4.0 * ue[t7 + 1][t6] + ue[t7 + 2][t6]);
                            }
                        }

                        for (int t6 = lbv; t6 <= ubv; t6++) {
                            forcing[t4][m_ - 3][t5][t6] = forcing[t4][m_ - 3][t5][t6] - dssp *
                                (ue[m_ - 3 - 2][t6] - 4.0 * ue[m_ - 3 - 1][t6] +
                                    6.0 * ue[m_ - 3][t6] - 4.0 * ue[m_ - 3 + 1][t6]);
                            
                            forcing[t4][m_ - 2][t5][t6] = forcing[t4][m_ - 2][t5][t6] - dssp *
                                (ue[m_ - 2 - 2][t6] - 4.0 * ue[m_ - 2 - 1][t6] + 5.0 * ue[m_ - 2][t6]);
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
