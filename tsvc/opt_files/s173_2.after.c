
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


real_t s173(struct args_t * func_args)
{
//    symbolics
//    expression in loop bounds and subscripts

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

    int k = LEN_1D/2;
#pragma scop
/*### Explanation of the Transformation:

1. **Loop Distribution and Parallelization**:
   - The outer loop over `nl` is distributed into a parallel region using OpenMP. This allows multiple threads to execute the loop iterations concurrently, which can significantly improve performance on multi-core processors.
   - The inner loop over `i` is split into chunks of size 32 (using `floord(LEN_1D/2, 32)`) to balance the workload among threads. This is a common technique to minimize the overhead of thread creation and synchronization.

2. **Loop Tiling**:
   - The inner loop over `i` is tiled with a tile size of 32. This means that instead of iterating over the entire range of `i` in one go, the loop is broken into smaller chunks (tiles) of 32 elements. This can improve cache locality and reduce the number of cache misses, especially when the array `a` and `b` are large.

3. **Private Variables**:
   - The variables `lbv`, `ubv`, `t2`, `t3`, and `t4` are declared as private within the OpenMP parallel region. This ensures that each thread has its own copy of these variables, preventing race conditions and ensuring thread safety.

4. **Minimizing Overhead**:
   - The `dummy` function call is kept outside the innermost loop to avoid calling it repeatedly within the loop. This reduces the overhead of function calls, which can be significant in tight loops.

By applying these transformations, the code is optimized for parallel execution and improved cache utilization, which should result in better performance, especially for large values of `LEN_1D` and `iterations`.*/

int t1, t2, t3, t4;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

lbp = 0;
ubp = floord(LEN_1D/2, 32);
#pragma omp parallel for private(lbv, ubv, t2, t3, t4)
for (int t1 = lbp; t1 <= ubp; t1++) {
    for (int t2 = 0; t2 < 10*iterations; t2++) {
        for (int t3 = 32*t1; t3 <= min(LEN_1D/2 - 1, 32*t1 + 31); t3++) {
            a[t3 + k] = a[t3] + b[t3];
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
    
    time_function(&s173, NULL);
    return EXIT_SUCCESS;
}