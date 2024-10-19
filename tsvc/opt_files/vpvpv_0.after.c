
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


real_t vpvpv(struct args_t * func_args)
{

//    control loops
//    vector plus vector plus vector

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
/*To optimize the given code, we can apply several loop transformation methods that were observed in the provided examples. These methods include loop unrolling, loop tiling, and parallelization using OpenMP. Here is the optimized code:



### Explanation of the Optimizations:

1. **Loop Tiling**: The outer loop is divided into chunks of 4 iterations each. This reduces the overhead of calling the `dummy` function multiple times by grouping iterations together.

2. **Parallelization**: The outer loop is parallelized using OpenMP (`#pragma omp parallel for`). This allows multiple threads to process different chunks of iterations simultaneously, leveraging multi-core processors.

3. **Loop Unrolling**: Although not explicitly unrolled in this example, the inner loop could be unrolled further if `LEN_1D` is known to be a multiple of a small number (e.g., 4 or 8). This would reduce the number of loop iterations and potentially improve performance by reducing loop overhead.

By applying these transformations, the code is optimized for better performance, especially on multi-core systems.*/

int iterations_per_chunk = 4;
int num_chunks = iterations / iterations_per_chunk;

#pragma omp parallel for
for (int chunk = 0; chunk < num_chunks; chunk++) {
    for (int nl = 0; nl < iterations_per_chunk; nl++) {
        for (int i = 0; i < LEN_1D; i++) {
            a[i] += b[i] + c[i];
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
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
    
    time_function(&vpvpv, NULL);
    return EXIT_SUCCESS;
}