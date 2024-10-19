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

#include <omp.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y) ((x) > (y)? (x) : (y))
#define min(x,y) ((x) < (y)? (x) : (y))

/* start param define */ 
#define N 100000
/* end parameters define */

/* start kernel func */
void Freebench_pifft1(double *out, int n)
{
    int j;
    int carry = 42;

    double time_start = omp_get_wtime();
#pragma scop
/**/

for (int iter = 0; iter < ITERATIONS; iter++){
    double temp = out[n];
    for (j = n; j > 2; j--) {
        out[j] = out[j - 1];
    }
    out[2] = carry;
}
#pragma endscop
    double time_end = omp_get_wtime();
    printf("%f\n", time_end - time_start);
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
