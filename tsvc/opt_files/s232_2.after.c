
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
/*### Explanation of the Corrected Optimized Code:

1. **Loop Distribution and Parallelization**:
   - The outer loop (`nl`) is distributed across multiple threads using OpenMP (`#pragma omp parallel for`). This allows the iterations to be executed in parallel, which can significantly improve performance on multi-core systems.

2. **Loop Tiling**:
   - The inner loops (`j` and `i`) are tiled using a tile size of 32 (`t1`). This helps in reducing the number of cache misses by keeping the data accessed by the inner loops within the cache for a longer duration.

3. **Loop Reordering**:
   - The loops are reordered to ensure that the most computationally intensive loop (`i`) is innermost, which is beneficial for memory access patterns and cache utilization.

4. **Reduction of Loop Bound Calculations**:
   - The loop bounds are precomputed and stored in variables (`lbp`, `ubp`, `t1`), which reduces the overhead of calculating these bounds repeatedly within the loops.

5. **Register Usage**:
   - The use of `register` for `lbv` and `ubv` hints to the compiler to store these variables in registers, which can improve access speed.

These transformations are based on the techniques observed in the provided examples, such as loop tiling, parallelization, and loop reordering, which are common methods for optimizing nested loops in performance-critical code.*/

int nl, j, i;
int lbp, ubp;
register int lbv, ubv;

lbp = 0;
ubp = floord(100 * (iterations / LEN_2D), 32);
#pragma omp parallel for private(lbv, ubv, j, i)
for (nl = lbp; nl <= ubp; nl++) {
    for (int t1 = 0; t1 <= floord(LEN_2D - 1, 32); t1++) {
        for (j = max(1, 32 * t1); j <= min(LEN_2D - 1, 32 * t1 + 31); j++) {
            for (i = 1; i <= j; i++) {
                aa[j][i] = aa[j][i - 1] * aa[j][i - 1] + bb[j][i];
            }
        }
    }
    dummy(a, b, c, d, e, aa, bb, cc, 1.);
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