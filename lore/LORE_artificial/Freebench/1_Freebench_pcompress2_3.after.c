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
void Freebench_pcompress2(double *freq, double *cum_freq, int No_of_symbols)
{
    int i;
    int cum = 0;

    double time_start = omp_get_wtime();
#pragma scop
   for (int iter = 0; iter < ITERATIONS; iter++){
       if (No_of_symbols > 0) {
           freq[No_of_symbols] = (freq[No_of_symbols] + 1) / 2;
           cum_freq[No_of_symbols] = cum;
           cum += freq[No_of_symbols];
       }
       for (i = No_of_symbols - 1; i > 0; i--) {
           freq[i] = (freq[i] + 1) / 2;
           cum_freq[i] = cum;
           cum += freq[i];
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

    Freebench_pcompress2(array_0, array_1, N);

    print_array_1d(array_0, N+1);

    free_array_1d(array_0, N+1);

    return 0;
}
