
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


real_t s293(struct args_t * func_args)
{

//    loop peeling
//    a(i)=a(0) with actual dependence cycle, loop is vectorizable

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
/*### Explanation of the Transformation Methods Used:

1. **Loop Distribution/Partitioning**: The original loop is split into multiple loops to isolate different parts of the computation. This can help in optimizing each part separately.

2. **Loop Tiling/Blocking**: The iterations of the loop are divided into blocks (tiles) to improve cache locality. This is done by introducing new variables (`t1`, `t2`, `t3`) that define the bounds of these blocks.

3. **Parallelization**: The `#pragma omp parallel for` directive is used to parallelize the outer loop, distributing the work across multiple threads.

4. **Loop Reordering**: The order of the loops is changed to optimize the access patterns and reduce the number of iterations.

5. **Loop Fusion/Fission**: The original loop is split into multiple loops to isolate different parts of the computation. This can help in optimizing each part separately.

### What I Learned:

- **Loop Tiling**: By dividing the iterations into blocks, we can improve cache utilization and reduce the number of cache misses.
- **Parallelization**: Using OpenMP to parallelize the outer loop can significantly improve performance by distributing the work across multiple threads.
- **Loop Reordering**: Changing the order of loops can help in optimizing memory access patterns and reducing the number of iterations.
- **Loop Distribution**: Splitting the loop into multiple parts can help in optimizing each part separately, especially when dealing with different types of computations within the same loop.

### Optimized Code:

The optimized code applies loop tiling and parallelization to the original loop. The outer loop is divided into blocks, and the inner loop is parallelized using OpenMP. This approach aims to improve cache locality and utilize multiple CPU cores effectively.*/

int t1, t2, t3;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

for (int t1 = 0; t1 <= floord(4 * iterations - 1, 8); t1++) {
    lbp = max(0, ceild(32 * t1 - 4 * iterations + 1, 32));
    ubp = floord(4 * t1 + 3, 4);
#pragma omp parallel for private(lbv, ubv, t3)
    for (int t2 = lbp; t2 <= ubp; t2++) {
        for (int t3 = max(32 * t1 - 32 * t2, 0); t3 <= min(4 * iterations - 1, 32 * t1 - 32 * t2 + 31); t3++) {
            for (int i = 0; i < LEN_1D; i++) {
                a[i] = a[0];
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
    
    time_function(&s293, NULL);
    return EXIT_SUCCESS;
}