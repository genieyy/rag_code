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

int t1, t2, t3, t4, t5, t6, t7, t8;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

for (int iter = 0; iter < ITERATIONS; iter++) {
    lbp = 0;
    ubp = floord(n_ - 1, 32);
    #pragma omp parallel for private(lbv, ubv, t3, t4, t5, t6, t7, t8)
    for (t2 = lbp; t2 <= ubp; t2++) {
        for (t3 = 0; t3 <= floord(q_ - 1, 32); t3++) {
            for (t4 = max(0, 32 * t2); t4 <= min(n_ - 1, 32 * t2 + 31); t4++) {
                lbv = 32 * t3;
                ubv = min(q_ - 1, 32 * t3 + 31);
                #pragma ivdep
                #pragma vector always
                for (t5 = lbv; t5 <= ubv; t5++) {
                    for (t6 = 0; t6 < p_; t6++) {
                        forcing[t4][1][t5][t6] = forcing[t4][1][t5][t6] - dssp *
                            (5.0 * ue[1][t6] - 4.0 * ue[1 + 1][t6] + ue[1 + 2][t6]);
                        
                        forcing[t4][2][t5][t6] = forcing[t4][2][t5][t6] - dssp *
                            (-4.0 * ue[2 - 1][t6] + 6.0 * ue[2][t6] -
                                4.0 * ue[2 + 1][t6] + ue[2 + 2][t6]);
                    }

                    for (t6 = 0; t6 < p_; t6++) {
                        for (t7 = 1 * 3; t7 <= m_ - 3 * 1 - 1; t7++) {
                            forcing[t4][t7][t5][t6] = forcing[t4][t7][t5][t6] - dssp *
                                (ue[t7 - 2][t6] - 4.0 * ue[t7 - 1][t6] +
                                    6.0 * ue[t7][t6] - 4.0 * ue[t7 + 1][t6] + ue[t7 + 2][t6]);
                        }
                    }
                    
                    for (t6 = 0; t6 < p_; t6++) {
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