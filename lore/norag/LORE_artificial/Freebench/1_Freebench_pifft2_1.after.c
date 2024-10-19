/*according to 
Benchmark: FreeBench default
Application: pifft
File: pifft.c
Function: mp_sprintf
Line: 1450
*/
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
#define N 200
#define M 300
/* end parameters define */

/* start kernel func */
void Freebench_pifft2(double *in, double *out, int n, int log10_radix)
{
    int x, y;

    double time_start = omp_get_wtime();
#pragma scop
/**/

for (int iter = 0; iter < ITERATIONS; iter++){
    for (int j = 3; j <= n + 1; j++) {
        x = in[j];
        double temp = x;
        for (int k = log10_radix - 1; k >= 0; k--) {
            y = (int)temp % 10;
            temp /= 10;
            out[(j - 3) * 3 + k] = (48 + y);
        }
    }
}
#pragma endscop
    double time_end = omp_get_wtime();
    printf("%f\n", time_end - time_start);
}
/* end kernel func */

int main(int argc, char *argv[])
{
    ARRAY_PREPARATION_1D(array_0, N+2);

    ARRAY_PREPARATION_1D(array_1, 3*N+M+1);

    Freebench_pifft2(array_0, array_1, N, M);

    print_array_1d(array_1, 3*N+M+1);

    free_array_1d(array_0, N+2);
    
    free_array_1d(array_1, 3*N+M+1);

    return 0;
}