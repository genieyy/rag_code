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
    const int TILE_SIZE = 8;
    for (int iter = 0; iter < ITERATIONS; iter++){
        for (int j_tile = 3; j_tile <= n + 1; j_tile += TILE_SIZE) {
            for (int j = j_tile; j < j_tile + TILE_SIZE && j <= n + 1; j++) {
                x = in[j];
                int base_index = (j - 3) * 3;
                for (int k = log10_radix - 1; k >= 0; k--) {
                    y = x % 10;
                    x /= 10;
                    out[base_index + k] = (48 + y);
                }
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