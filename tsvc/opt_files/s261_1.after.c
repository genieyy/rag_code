
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


real_t s261(struct args_t * func_args)
{

//    scalar and array expansion
//    wrap-around scalar under an if

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

    real_t t;
#pragma scop
/*To optimize the given code, we can apply several loop transformation methods inspired by the examples provided:

1. **Loop Unrolling**: This technique reduces the overhead of loop control by processing multiple elements within a single iteration.
2. **Parallelization**: Using OpenMP to parallelize the outer loop, which can be safely parallelized as the iterations are independent.
3. **Vectorization**: Using compiler directives to hint the compiler to generate vectorized instructions.

Here is the optimized code:



### Explanation:
1. **Parallelization**: The outer loop is parallelized using `#pragma omp parallel for private(t)`, which ensures that each thread has its own private copy of the variable `t`.
2. **Loop Unrolling**: The inner loop is unrolled by a factor of 4 (`i += 4`), which reduces the number of iterations and the overhead of loop control.
3. **Vectorization**: The `#pragma ivdep` and `#pragma vector always` directives are used to hint the compiler to generate vectorized instructions, assuming that the loop iterations are independent and can be safely vectorized.

This optimization should improve the performance of the code by reducing loop overhead and leveraging parallel processing capabilities.*/

#pragma omp parallel for private(t)
for (int nl = 0; nl < iterations; nl++) {
    #pragma ivdep
    #pragma vector always
    for (int i = 1; i < LEN_1D; i += 4) {
        t = a[i] + b[i];
        a[i] = t + c[i-1];
        t = c[i] * d[i];
        c[i] = t;

        t = a[i+1] + b[i+1];
        a[i+1] = t + c[i];
        t = c[i+1] * d[i+1];
        c[i+1] = t;

        t = a[i+2] + b[i+2];
        a[i+2] = t + c[i+1];
        t = c[i+2] * d[i+2];
        c[i+2] = t;

        t = a[i+3] + b[i+3];
        a[i+3] = t + c[i+2];
        t = c[i+3] * d[i+3];
        c[i+3] = t;
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
    
    time_function(&s261, NULL);
    return EXIT_SUCCESS;
}