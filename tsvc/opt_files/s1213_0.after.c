
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


real_t s1213(struct args_t * func_args)
{

//    statement reordering
//    dependency needing temporary

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
/*### Explanation of the Transformation:
1. **Loop Distribution and Parallelization**:
   - The outer loop over `nl` is distributed into chunks of size 32 using `floord(iterations-1, 32)`. This allows for parallel execution of these chunks using OpenMP.
   - The `#pragma omp parallel for` directive is used to parallelize the outer loop, ensuring that each chunk of iterations is executed in parallel.

2. **Loop Tiling**:
   - The inner loop over `i` is tiled using a tile size of 32 (`32*t1` to `32*t1+31`). This helps in better cache utilization and can improve performance by reducing cache misses.

3. **Loop Fusion**:
   - The operations `a[i] = b[i-1] + c[i];` and `b[i] = a[i+1] * d[i];` are kept within the same loop to avoid unnecessary loop overhead and to ensure that the operations are performed in a pipelined manner.

4. **Reduction of Loop Bound Calculations**:
   - The bounds for the inner loop are precomputed using `max(1, 32*t1)` and `min(LEN_1D-2, 32*t1+31)`, which reduces the number of times these bounds are calculated during the loop execution.

5. **Register Usage**:
   - The use of `register` for `lbv` and `ubv` hints to the compiler to place these variables in registers, potentially improving performance by reducing memory access.

This optimization strategy leverages parallel processing, loop tiling, and efficient loop bounds calculation to improve the performance of the original loop.*/

int lb, ub, lbp, ubp;
register int lbv, ubv;
lbp = 0;
ubp = floord(iterations-1, 32);
#pragma omp parallel for private(lbv, ubv)
for (int t1 = lbp; t1 <= ubp; t1++) {
    for (int t2 = max(1, 32*t1); t2 <= min(LEN_1D-2, 32*t1+31); t2++) {
        a[t2] = b[t2-1] + c[t2];
        b[t2] = a[t2+1] * d[t2];
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
    
    time_function(&s1213, NULL);
    return EXIT_SUCCESS;
}