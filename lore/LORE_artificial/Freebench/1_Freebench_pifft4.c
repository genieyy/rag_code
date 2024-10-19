/*according to 
Benchmark: FreeBench default
Application: pifft
File: pifft.c
Function: mp_unsgn_imul
Line: 662
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
void Freebench_pifft4(double *out, int n)
{
    int shift = 5;

#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++){
    for (int j = n; j >= shift + 1; j--) {
        out[j] = out[j - shift];
    }
}
#pragma endscop
}
/* end kernel func */

int main(int argc, char *argv[])
{
    ARRAY_PREPARATION_1D(array_0, N+1);

    Freebench_pifft4(array_0, N);

    print_array_1d(array_0, N+1);

    free_array_1d(array_0, N+1);

    return 0;
}