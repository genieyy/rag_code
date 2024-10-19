/*according to 
Benchmark: FreeBench default
Application: pifft
File: pifft.c
Function: mp_unexp_add
Line: 500
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
#define N 20000
/* end parameters define */

/* start kernel func */
void Freebench_pifft3(double *in1, double *in2, double *out, int n)
{
    int x;
    int carry = 0;
    int radix = 10;

    double time_start = omp_get_wtime();
#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++){
    for (int j = n - 1; j > 0; j -= 4) {
        double x1 = in1[j - 1] + in2[j - 1] - carry;
        carry = (x1 >= radix) ? -1 : 0;
        out[j] = x1 - (radix & carry);
        
        double x2 = in1[j - 2] + in2[j - 2] - carry;
        carry = (x2 >= radix) ? -1 : 0;
        out[j - 1] = x2 - (radix & carry);
        
        double x3 = in1[j - 3] + in2[j - 3] - carry;
        carry = (x3 >= radix) ? -1 : 0;
        out[j - 2] = x3 - (radix & carry);
        
        double x4 = in1[j - 4] + in2[j - 4] - carry;
        carry = (x4 >= radix) ? -1 : 0;
        out[j - 3] = x4 - (radix & carry);
    }
}
#pragma endscop
    double time_end = omp_get_wtime();
    printf("%f\n", time_end - time_start);
}
/* end kernel func */

int main(int argc, char *argv[])
{
    ARRAY_PREPARATION_1D(array_0, N+1);
    ARRAY_PREPARATION_1D(array_1, N+1);
    ARRAY_PREPARATION_1D(array_2, N+1);

    Freebench_pifft3(array_0, array_1, array_2, N);

    print_array_1d(array_2, N+1);

    free_array_1d(array_0, N+1);
    free_array_1d(array_1, N+1);
    free_array_1d(array_2, N+1);

    return 0;
}