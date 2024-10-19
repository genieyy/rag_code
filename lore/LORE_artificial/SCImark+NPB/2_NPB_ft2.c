#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lore.h"

/* start param define */
#define N 100000
/* end parameters define */

/* start kernel func */
void NPB_ft2(double *ex, int EXPMAX)
{
    int i;

#pragma scop
for (int iter = 0; iter < ITERATIONS; iter++){
    for (i = 2; i <= EXPMAX; i++)
    {
        ex[i] = ex[i - 1] * ex[1];
    }
}
#pragma endscop
}
/* end kernel func */

int main(int argc, char *argv[])
{
    ARRAY_PREPARATION_1D(array_0, N + 1);

    NPB_ft2(array_0, N);

    print_array_1d(array_0, N + 1);

    free_array_1d(array_0, N + 1);

    return 0;
}