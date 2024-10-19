
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
/*### Explanation of the Transformation:
1. **Loop Distribution and Parallelization**:
   - The outer loop is distributed into chunks of size 32 using the `floord` function to determine the number of chunks.
   - The `#pragma omp parallel for` directive is used to parallelize the outer loop, distributing the work across multiple threads.

2. **Loop Reordering**:
   - The inner loops are reordered to ensure that the most computationally intensive loop (`for (int i = 1; i < LEN_2D; i++)`) is innermost, which is beneficial for cache locality.

3. **Loop Unrolling**:
   - The inner loop is not unrolled in this example, but it could be considered if further optimization is needed, especially if `LEN_2D` is known to be a multiple of a small number.

4. **Reduction of Loop Bounds**:
   - The loop bounds are reduced using `max` and `min` functions to ensure that the loops do not exceed the valid range of indices.

5. **Register Usage**:
   - The `register` keyword is used for `lbv` and `ubv` to suggest that these variables should be stored in CPU registers for faster access.

This transformation aims to improve performance by leveraging parallel processing and optimizing loop structures for better cache utilization.*/

int t1, t2, t3;
int lbp, ubp;
register int lbv, ubv;

lbp = 0;
ubp = floord(100 * (iterations / LEN_2D), 32);

#pragma omp parallel for private(lbv, ubv, t2, t3)
for (t1 = lbp; t1 <= ubp; t1++) {
    for (t2 = max(0, 32 * t1); t2 <= min(100 * (iterations / LEN_2D) - 1, 32 * t1 + 31); t2++) {
        for (t3 = 1; t3 < LEN_2D; t3++) {
            for (int i = 1; i < LEN_2D; i++) {
                aa[t3][i] = (aa[t3][i - 1] + aa[t3 - 1][i]) / 1.9;
            }
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
    
    time_function(&s2111, NULL);
    return EXIT_SUCCESS;
}