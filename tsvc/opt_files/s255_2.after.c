
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


real_t s255(struct args_t * func_args)
{

//    scalar and array expansion
//    carry around variables, 2 levels

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

    real_t x, y;
#pragma scop
/*### Explanation of Optimizations:
1. **Reduction of Array Accesses**: By storing `b[LEN_1D-1]` and `b[LEN_1D-2]` in variables `x` and `y` before the loop, we reduce the number of array accesses. This is particularly beneficial if array accesses are expensive (e.g., due to cache misses).

2. **Loop Unrolling**: Although not explicitly unrolled, the loop structure is kept simple to allow for potential compiler optimizations. The loop is already quite simple, so further unrolling might not provide significant benefits without additional context.

3. **Constant Factor Extraction**: The constant multiplication factor `(real_t).333` is extracted outside the loop and stored in a variable `factor`. This avoids recalculating the constant during each iteration, which can be slightly more efficient.

4. **Temporary Variable Usage**: A temporary variable `temp` is used to store `b[i]` before it is assigned to `x`. This avoids the need to access `b[i]` multiple times within the loop, which can be beneficial for performance.

These optimizations are based on the principles of reducing redundant computations and minimizing memory accesses, which are common techniques in loop optimization.*/

for (int nl = 0; nl < iterations; nl++) {
    double x = b[LEN_1D-1];
    double y = b[LEN_1D-2];
    double factor = (real_t).333;
    for (int i = 0; i < LEN_1D; i++) {
        double temp = b[i];
        a[i] = (temp + x + y) * factor;
        y = x;
        x = temp;
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
    
    time_function(&s255, NULL);
    return EXIT_SUCCESS;
}