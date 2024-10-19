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

/* start param define */ 
#define N 200
#define M 300
/* end parameters define */

/* start kernel func */
void Freebench_pifft2(double *in, double *out, int n, int log10_radix)
{
    int x, y;

#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++){
    for (int j = 3; j <= n + 1; j++) {
        x = in[j];
        for (int k = log10_radix - 1; k >= 0; k--) {
            y = x % 10;
            x /= 10;
            out[(j - 3) * 3 + k] = (48 + y);
        }
    }
}
#pragma endscop
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