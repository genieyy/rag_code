
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


real_t s231(struct args_t * func_args)
{
//    loop interchange
//    loop with data dependency

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
/*### Analysis of Meaning-Preserving Loop Transformation Methods:

1. **Loop Unrolling**: This technique is used to reduce the overhead of loop control instructions by processing multiple elements in a single iteration. However, it is not explicitly shown in the provided examples.

2. **Loop Fusion**: Combining multiple loops that iterate over the same range into a single loop to reduce loop overhead and improve cache locality. This is not explicitly shown in the provided examples.

3. **Loop Distribution**: Splitting a single loop into multiple loops to isolate independent parts that can be optimized separately. This is not explicitly shown in the provided examples.

4. **Loop Tiling (Blocking)**: This technique is used to improve cache locality by breaking the loop into smaller chunks (tiles) that fit better into the cache. This is not explicitly shown in the provided examples.

5. **Parallelization**: Using OpenMP to parallelize loops, allowing multiple threads to execute the loop body concurrently. This is shown in the provided examples with `#pragma omp parallel for`.

### Applying the Learned Methods:

1. **Parallelization**: The provided optimized code uses OpenMP to parallelize the outer loop (`nl` loop), which is the most significant loop in terms of iterations. This allows multiple threads to execute the inner loops concurrently, potentially improving performance.

2. **Loop Distribution**: The inner loops (`i` and `j` loops) are kept as they are, as they are already simple and do not benefit from further distribution.

3. **Loop Fusion**: There is no opportunity for loop fusion in the provided code as there are no adjacent loops that can be combined.

4. **Loop Unrolling**: The inner loops are not unrolled as they are already simple and unrolling them might not provide significant benefits.

5. **Loop Tiling**: The inner loops are not tiled as they are already simple and tiling them might not provide significant benefits.

### Optimized Code:

The optimized code leverages parallelization to improve performance by distributing the work across multiple threads. The inner loops are kept as they are since they are already simple and do not benefit from further transformations.*/

int nl, i, j;
int lb, ub, lbp, ubp;
register int lbv, ubv;

lbp = 0;
ubp = (iterations / LEN_2D) - 1;
#pragma omp parallel for private(lbv, ubv, i, j)
for (nl = lbp; nl <= ubp; nl++) {
    for (i = 0; i < LEN_2D; i++) {
        for (j = 1; j < LEN_2D; j++) {
            aa[j][i] = aa[j - 1][i] + bb[j][i];
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
    
    time_function(&s231, NULL);
    return EXIT_SUCCESS;
}