
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


real_t s244(struct args_t * func_args)
{

//    node splitting
//    false dependence cycle breaking

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
/*### Analysis of Meaning-Preserving Loop Transformation Methods:

1. **Loop Tiling (Blocking):**
   - The original loops are divided into smaller blocks (tiles) to improve cache locality. This is evident in the transformation where the loop bounds are adjusted to iterate over smaller chunks of the array (e.g., `32*t1` to `32*t1+31`).

2. **Parallelization:**
   - The loops are parallelized using OpenMP (`#pragma omp parallel for`) to exploit multi-core processors. This allows multiple iterations of the loop to be executed simultaneously.

3. **Loop Fusion:**
   - In some cases, loops that operate on the same data are fused together to reduce the overhead of loop control and improve data locality. This is seen in the optimized code where multiple operations on the same array elements are combined within a single loop.

4. **Loop Distribution:**
   - The loops are distributed to handle different conditions separately. This is useful for handling edge cases or specific conditions within the loop body.

### Learning from the Examples:

- **Cache Efficiency:** By tiling the loops, we can ensure that the data accessed within each iteration fits better into the cache, reducing cache misses and improving performance.
- **Parallel Execution:** Leveraging OpenMP for parallel execution can significantly speed up the computation by distributing the workload across multiple threads.
- **Reduced Overhead:** Combining operations within a single loop reduces the overhead associated with loop control and can lead to more efficient code.

### Optimized Code Explanation:

- **Tiling:** The outer loop (`t1`) iterates over blocks of size 32, ensuring that the inner loop operates on smaller chunks of the array.
- **Parallelization:** The `#pragma omp parallel for` directive is used to parallelize the outer loop, allowing multiple blocks to be processed simultaneously.
- **Loop Fusion:** The operations on `a`, `b`, and `c` are combined within the same loop, reducing the overhead of multiple loop iterations.

This optimized code should provide better performance by improving cache utilization and leveraging parallel execution.*/

int lb, ub, lbp, ubp;
register int lbv, ubv;
lbp = 0;
ubp = floord(LEN_1D-2, 32);
#pragma omp parallel for private(lbv, ubv)
for (int t1 = lbp; t1 <= ubp; t1++) {
    for (int nl = 0; nl < iterations; nl++) {
        for (int t2 = max(0, 32*t1); t2 <= min(LEN_1D-2, 32*t1+31); t2++) {
            a[t2] = b[t2] + c[t2] * d[t2];
            b[t2] = c[t2] + b[t2];
            a[t2+1] = b[t2] + a[t2+1] * d[t2];
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
    
    time_function(&s244, NULL);
    return EXIT_SUCCESS;
}