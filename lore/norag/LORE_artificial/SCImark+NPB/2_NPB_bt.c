#include <stdio.h>  
#include <stdlib.h>  
#include <string.h>
#include <math.h>
#include "lore.h"

/* start param define */ 
#define N 20
#define M 30
#define Q 40
#define P 6
/* end parameters define */

/* start kernel func */
void NPB_bt(double ****u, double ****rhs, int n, int m, int q, int p) {
    int i, j, k, l;
#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++){
    for (i = 1; i < n - 1 ; i++) {
        for (j = 1; j < m - 1; j++) {
            for (k = 1; k < q - 1; k++) {
                for (l = 0; l < p - 1; l++) {
                    u[i][j][k][l] = u[i][j][k][l] + rhs[i][j][k][l];
                }
            }
        }
    }
}
#pragma endscop
}
/* end kernel func */

int main() {
    ARRAY_PREPARATION_4D(array_0, N + 1, M + 1, Q + 1, P + 1);
    ARRAY_PREPARATION_4D(array_1, N + 1, M + 1, Q + 1, P + 1);

    NPB_bt(array_0, array_1, N, M, Q, P);
    
    print_array_4d(array_0, N, M, Q, P);

    free_array_4d(array_0, N + 1, M + 1, Q + 1, P + 1);
    free_array_4d(array_1, N + 1, M + 1, Q + 1, P + 1);

    return 0;
}