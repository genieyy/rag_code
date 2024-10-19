
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


real_t s113(struct args_t * func_args)
{

//    linear dependence testing
//    a(i)=a(1) but no actual dependence cycle

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
/*### Explanation of the Transformation:

1. **Loop Distribution and Parallelization**:
   - The outer loop over `nl` remains unchanged as it iterates over the number of iterations.
   - The inner loop over `i` is parallelized using OpenMP to distribute the workload across multiple threads.
   - The loop bounds are adjusted to ensure that each thread processes a contiguous block of indices, which can improve cache locality and reduce false sharing.

2. **Loop Tiling**:
   - The inner loop over `i` is tiled with a tile size of 32, which is a common choice for cache optimization. This ensures that each thread processes a block of 32 elements at a time, which can help in better utilization of the CPU cache.

3. **Private Variables**:
   - The variables `lbv`, `ubv`, and `i` are declared as private within the OpenMP parallel region to ensure that each thread has its own copy of these variables, preventing race conditions.

4. **Loop Bounds Adjustment**:
   - The loop bounds are adjusted using `max` and `min` functions to ensure that the loop iterates over valid indices within the block assigned to each thread.

This transformation aims to improve the performance of the original loop by leveraging parallel execution and optimizing cache usage through loop tiling.*/

int nl, i;
int lb, ub, lbp, ubp;
register int lbv, ubv;

for (nl = 0; nl < 4 * iterations; nl++) {
    lbp = 0;
    ubp = floord(LEN_1D - 1, 32);
#pragma omp parallel for private(lbv, ubv, i)
    for (lbv = lbp; lbv <= ubp; lbv++) {
        for (i = max(1, 32 * lbv); i <= min(LEN_1D - 1, 32 * lbv + 31); i++) {
            a[i] = a[0] + b[i];
        }
    }
    dummy(a, b, c, d, e, aa, bb, cc, 0.);
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
    
    time_function(&s113, NULL);
    return EXIT_SUCCESS;
}