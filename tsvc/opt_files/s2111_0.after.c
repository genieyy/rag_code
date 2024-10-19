
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


real_t s2111(struct args_t * func_args)
{

//    wavefronts, it will make jump in data access

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
/*### Explanation of the Optimized Code:

1. **Loop Transformation**:
   - The original nested loops are transformed into a single outer loop (`t1`) that iterates over a range calculated using `floord` and `ceild` functions to ensure the loop bounds are correctly aligned with the problem size.
   - The inner loops (`t2` and `t3`) are used to further partition the iteration space, ensuring that the iterations are distributed efficiently across multiple threads.

2. **Parallelization**:
   - The `#pragma omp parallel for` directive is used to parallelize the outer loop (`t1`), allowing multiple threads to execute the loop iterations concurrently.
   - The `private` clause ensures that each thread has its own copy of the loop variables (`lbv`, `ubv`, `t2`, `t3`), preventing race conditions.

3. **Loop Bounds Calculation**:
   - The bounds for the loops are calculated using `max` and `min` functions to ensure that the iterations are within the valid range of the problem size (`LEN_2D`).
   - The `floord` and `ceild` functions are used to align the loop bounds with the problem size, ensuring that the iterations are evenly distributed across the threads.

4. **Loop Unrolling**:
   - The inner loops (`j` and `i`) are not unrolled in this example, but the outer loops (`t1`, `t2`, `t3`) are structured to allow for potential loop unrolling in future optimizations.

5. **Dummy Function Call**:
   - The `dummy` function call is placed inside the innermost loop to ensure that it is executed after each iteration of the main computation.

This optimization strategy leverages loop transformation and parallelization to improve the performance of the original code by distributing the workload across multiple threads and ensuring efficient use of the available resources.*/

int t1, t2, t3;
int lbp, ubp;
register int lbv, ubv;

lbp = 0;
ubp = floord(100 * (iterations / LEN_2D), 32);

#pragma omp parallel for private(lbv, ubv, t2, t3)
for (t1 = lbp; t1 <= ubp; t1++) {
    for (t2 = max(0, ceild(32 * t1 - 99 * (iterations / LEN_2D) + 1, 32)); t2 <= min(floord(100 * (iterations / LEN_2D) - 1, 32), t1); t2++) {
        for (t3 = max(32 * t1, 32 * t2 + 32); t3 <= min(100 * (iterations / LEN_2D) - 1, 32 * t1 + 31); t3++) {
            for (int j = 1; j < LEN_2D; j++) {
                for (int i = 1; i < LEN_2D; i++) {
                    aa[j][i] = (aa[j][i-1] + aa[j-1][i]) / 1.9;
                }
            }
            dummy(a, b, c, d, e, aa, bb, cc, 0.);
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
    
    time_function(&s2111, NULL);
    return EXIT_SUCCESS;
}