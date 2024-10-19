/*according to 
Benchmark: FreeBench default
Application: pifft
File: pifft.c
Function: mp_mul_d2i
Line: 1053
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lore.h"

/* start param define */ 
#define N 100000
/* end parameters define */

/* start kernel func */
void Freebench_pifft1(double *out, int n)
{
    int j;
    int carry = 42;

#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++){
    for (j = n + 1; j > 2; j--) {
        out[j] = out[j - 1];
    }
    out[2] = carry;
}
#pragma endscop
}
/* end kernel func */

int main(int argc, char *argv[])
{
    ARRAY_PREPARATION_1D(array_0, N+2);

    Freebench_pifft1(array_0, N);

    print_array_1d(array_0, N+2);

    free_array_1d(array_0, N+2);

    return 0;
}
