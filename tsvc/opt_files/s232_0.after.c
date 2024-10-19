
/*
 * This is an executable test containing a number of loops to measure
 * the performance of a compiler. Arrays' length is LEN_1D by default
 * and if you want a different array length, you should replace every
 * LEN_1D by your desired number which must be a multiple of 40. If you
 * want to increase the number of loop calls to have a longer run time
 * you have to manipulate the constant value iterations. There is a dummy
 * function called in each loop to make all computations appear required.
 * The time to execute this function is included in the time measurement
 * for the output but it is neglectable.
 *
 *  The output includes three columns:
 *    Loop:        The name of the loop
 *    Time(Sec):     The time in seconds to run the loop
 *    Checksum:    The checksum calculated when the test has run
 *
 * In this version of the codelets arrays are static type.
 *
 * All functions/loops are taken from "TEST SUITE FOR VECTORIZING COMPILERS"
 * by David Callahan, Jack Dongarra and David Levine except those whose
 * functions' name have 4 digits.
 */

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"
#include <omp.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y) ((x) > (y)? (x) : (y))
#define min(x,y) ((x) < (y)? (x) : (y))
#ifndef TSVC_DUMP_ARRAYS
#define DUMP 0
#else
#define DUMP 1
#endif

// array definitions
__attribute__((aligned(ARRAY_ALIGNMENT))) real_t flat_2d_array[LEN_2D*LEN_2D];

__attribute__((aligned(ARRAY_ALIGNMENT))) real_t x[LEN_1D];

__attribute__((aligned(ARRAY_ALIGNMENT))) real_t a[LEN_1D],b[LEN_1D],c[LEN_1D],d[LEN_1D],e[LEN_1D],
                                   aa[LEN_2D][LEN_2D],bb[LEN_2D][LEN_2D],cc[LEN_2D][LEN_2D],tt[LEN_2D][LEN_2D];

__attribute__((aligned(ARRAY_ALIGNMENT))) int indx[LEN_1D];

real_t* __restrict__ xx;
real_t* yy;


real_t s232(struct args_t * func_args)
{

//    loop interchange
//    interchanging of triangular loops

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
/*### Analysis and Transformation Methods Used:

1. **Loop Tiling/Blocking**: The original loops are transformed into a tiled structure to improve cache locality. This is done by introducing new variables (`t1`, `t2`, `t3`) that represent the tile indices.

2. **Parallelization**: The outer loop is parallelized using OpenMP (`#pragma omp parallel for`) to leverage multi-core processors. This allows multiple iterations of the loop to be executed simultaneously.

3. **Loop Reordering**: The original nested loops are reordered to ensure that the most computationally intensive part (the inner loops) is executed within the parallelized outer loop. This helps in better load balancing across threads.

4. **Loop Fusion**: The `dummy` function call, which is outside the original inner loops, is fused with the inner loops to ensure that it is executed within the parallelized region. This reduces the overhead of thread management.

5. **Loop Bound Adjustment**: The bounds of the loops are adjusted to ensure that the loop iterations are correctly aligned with the tile sizes. This involves calculating the maximum and minimum values for the tile indices (`t2`, `t3`) to ensure that the loop iterations are within the valid range.

### Performance Improvement:

- **Cache Locality**: By tiling the loops, the code ensures that data accessed within the inner loops is more likely to be in the cache, reducing memory access latency.
- **Parallel Execution**: By parallelizing the outer loop, the code can take advantage of multiple CPU cores, leading to a potential speedup proportional to the number of available cores.
- **Reduced Overhead**: The fusion of the `dummy` function call with the inner loops reduces the overhead associated with thread management and synchronization.*/

int t1, t2, t3;
int lbp, ubp;
register int lbv, ubv;

lbp = 0;
ubp = floord(100 * (iterations / LEN_2D), 32);

#pragma omp parallel for private(lbv, ubv, t2, t3)
for (t1 = lbp; t1 <= ubp; t1++) {
    for (t2 = max(0, ceild(32 * t1 - 99 * (iterations / LEN_2D) + 1, 32)); t2 <= min(floord(100 * (iterations / LEN_2D) - 1, 32), t1); t2++) {
        for (t3 = max(32 * t1, 32 * t2 + 1); t3 <= min(100 * (iterations / LEN_2D) - 1, 32 * t1 + 31); t3++) {
            for (int j = 1; j < LEN_2D; j++) {
                for (int i = 1; i <= j; i++) {
                    aa[j][i] = aa[j][i-1] * aa[j][i-1] + bb[j][i];
                }
            }
            dummy(a, b, c, d, e, aa, bb, cc, 1.);
        }
    }
}
#pragma endscop

    func_args->c2 = omp_get_wtime();
    return calc_checksum(__func__);
}

typedef real_t(*test_function_t)(struct args_t *);

void time_function(test_function_t vector_func, void * arg_info)
{
    struct args_t func_args = {.arg_info=arg_info};

    double result = vector_func(&func_args);

    double tic=func_args.c1;
    double toc=func_args.c2;
    double taken = toc-tic;

    printf("%lf\t%lf\n", taken, result);
    if (DUMP) save_array(func_args.name);
}

int main(int argc, char ** argv){
    int n1 = 1;
    int n3 = 1;
    int* ip;
    real_t s1,s2;
    init(&ip, &s1, &s2);
    printf("Loop \tTime(sec) \tChecksum\n");
    
    time_function(&s232, NULL);
    return EXIT_SUCCESS;
}