
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


real_t s281(struct args_t * func_args)
{

//    crossing thresholds
//    index set splitting
//    reverse data access

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

    real_t x;
#pragma scop
/*To optimize the given code, we can apply several loop transformation methods inspired by the examples provided:

1. **Loop Interchange**: This can be used to change the order of nested loops to improve cache locality.
2. **Loop Unrolling**: This can be used to reduce loop overhead and improve instruction-level parallelism.
3. **Vectorization**: This can be used to exploit SIMD (Single Instruction, Multiple Data) instructions for better performance.

Here is the optimized code:



### Explanation:
1. **Loop Unrolling**: The inner loop is unrolled by a factor of 4 to reduce the number of loop iterations and improve instruction-level parallelism.
2. **Register Usage**: Values are loaded into registers (`a_val1`, `b_val1`, `c_val1`, etc.) to avoid repeated memory accesses, which can be costly.
3. **Vectorization**: Although not explicitly using SIMD instructions, the unrolled loop structure can be beneficial for compilers that can auto-vectorize the code.

This optimization should improve the performance of the code by reducing loop overhead and improving cache locality.*/

for (int nl = 0; nl < iterations; nl++) {
    for (int i = 0; i < LEN_1D; i += 4) {
        // Load values into registers
        real_t a_val1 = a[LEN_1D-i-1];
        real_t a_val2 = a[LEN_1D-i-2];
        real_t a_val3 = a[LEN_1D-i-3];
        real_t a_val4 = a[LEN_1D-i-4];

        real_t b_val1 = b[i];
        real_t b_val2 = b[i+1];
        real_t b_val3 = b[i+2];
        real_t b_val4 = b[i+3];

        real_t c_val1 = c[i];
        real_t c_val2 = c[i+1];
        real_t c_val3 = c[i+2];
        real_t c_val4 = c[i+3];

        // Compute x values
        real_t x1 = a_val1 + b_val1 * c_val1;
        real_t x2 = a_val2 + b_val2 * c_val2;
        real_t x3 = a_val3 + b_val3 * c_val3;
        real_t x4 = a_val4 + b_val4 * c_val4;

        // Update arrays
        a[i] = x1 - (real_t)1.0;
        a[i+1] = x2 - (real_t)1.0;
        a[i+2] = x3 - (real_t)1.0;
        a[i+3] = x4 - (real_t)1.0;

        b[i] = x1;
        b[i+1] = x2;
        b[i+2] = x3;
        b[i+3] = x4;
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
    
    time_function(&s281, NULL);
    return EXIT_SUCCESS;
}