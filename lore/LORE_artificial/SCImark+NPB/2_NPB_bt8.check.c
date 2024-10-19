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
#define N 40
#define M 50
#define Q 60
#define P 5
/* end parameters define */

/* start kernel func */
void NPB_bt8(double ****forcing, int n_, int m_, int q_, int p_) {
    int i, j, k, m;
    double dssp = 0.5; // Example value, replace with actual value

    double time_start = omp_get_wtime();
 #pragma scop
for (int iter = 0; iter < ITERATIONS; iter++){
    for (i = 1; i < n_ - 1; i++) {
        for (j = 1; j < m_ - 1; j++) {
            for (k = 1; k < q_ - 1; k++) {
                for (m = 0; m < p_; m++) {
                    forcing[i][j][k][m] = -1.0 * forcing[i][j][k][m];
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
    ARRAY_PREPARATION_4D(array_0, N+1, M+1, Q+1, P+1);

    NPB_bt8(array_0, N, M, Q, P);

    print_array_4d(array_0, N+1, M+1, Q+1, P+1);

    free_array_4d(array_0, N+1, M+1, Q+1, P+1);
    
    return 0;
}