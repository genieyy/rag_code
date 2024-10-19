
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


real_t s114(struct args_t * func_args)
{

//    linear dependence testing
//    transpose vectorization
//    Jump in data access - not vectorizable

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
/*### Analysis of Meaning-Preserving Loop Transformation Methods

1. **Loop Unrolling**: The original code has nested loops, and the optimized code does not explicitly unroll loops. However, the use of `#pragma omp parallel for` can implicitly lead to loop unrolling by the compiler, especially when combined with OpenMP's parallelization.

2. **Parallelization**: The optimized code uses OpenMP to parallelize the outer loop (`nl` loop). This allows multiple threads to execute the loop iterations concurrently, which can significantly improve performance on multi-core systems.

3. **Loop Fusion**: The optimized code does not explicitly fuse loops. However, the parallelization of the outer loop can be seen as a form of loop fusion where the outer loop is split into multiple threads, each handling a subset of the iterations.

4. **Loop Distribution**: The optimized code does not explicitly distribute loops. However, the use of OpenMP implicitly distributes the loop iterations across multiple threads.

5. **Loop Tiling**: The optimized code does not explicitly tile loops. However, the use of OpenMP can lead to a form of implicit tiling where different parts of the loop are executed by different threads.

### Learning from the Examples

- **Parallelization**: The examples show that using OpenMP can significantly improve performance by parallelizing loops. This is particularly useful for loops with large iteration counts.
- **Loop Nesting**: The examples demonstrate that careful nesting of loops can help in optimizing memory access patterns and reducing cache misses.
- **Register Usage**: The use of `register` variables in the examples suggests that using register variables for frequently accessed variables can improve performance by reducing memory access overhead.

### Optimized Code Explanation

- **Parallelization**: The outer loop (`nl` loop) is parallelized using OpenMP, allowing multiple threads to execute the loop iterations concurrently.
- **Loop Nesting**: The inner loops (`i` and `j` loops) remain nested as in the original code, ensuring that the memory access patterns are preserved.
- **Register Usage**: Although not explicitly shown in the optimized code, the use of `register` variables can be considered for frequently accessed variables to improve performance.

This optimized code leverages parallelization to improve performance while preserving the original meaning of the loop structure.*/

int nl, i, j;
int lbp, ubp;
register int lbv, ubv;

lbp = 0;
ubp = (iterations / LEN_2D) - 1;

#pragma omp parallel for private(lbv, ubv, i, j)
for (nl = lbp; nl <= ubp; nl++) {
    for (i = 0; i < LEN_2D; i++) {
        for (j = 0; j < i; j++) {
            aa[i][j] = aa[j][i] + bb[i][j];
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
    
    time_function(&s114, NULL);
    return EXIT_SUCCESS;
}