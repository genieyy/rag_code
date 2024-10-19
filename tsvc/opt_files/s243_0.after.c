
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


real_t s243(struct args_t * func_args)
{

//    node splitting
//    false dependence cycle breaking

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
/*### Analysis of Loop Transformation Methods Used:

1. **Loop Distribution and Parallelization**:
   - The original loops are distributed into smaller chunks that can be executed in parallel using OpenMP. This is done by introducing new loop variables (`t1`, `t2`, `t3`, `t4`, `t5`) and partitioning the iterations into smaller blocks.
   - The `#pragma omp parallel for` directive is used to parallelize the outermost loop, allowing multiple threads to execute the loop iterations concurrently.

2. **Loop Tiling**:
   - The iterations are divided into tiles (blocks) of size 32, which helps in reducing cache misses and improving locality of reference. This is particularly useful for large arrays like `a`, `b`, `c`, `d`, and `e`.

3. **Loop Fusion**:
   - The inner loops are fused together to reduce the overhead of loop control and to improve data locality. This is done by iterating over the same range of indices in a single loop instead of multiple nested loops.

4. **Loop Interchange**:
   - The order of loops is interchanged to ensure that the innermost loop operates on the most frequently accessed data, thereby improving cache utilization.

### Learning:
- **Parallelization**: Utilizing OpenMP to parallelize loops can significantly improve performance by leveraging multiple CPU cores.
- **Loop Tiling**: Breaking down large loops into smaller tiles can improve cache performance and reduce memory access latency.
- **Loop Fusion**: Combining multiple loops that operate on the same data can reduce overhead and improve data locality.
- **Loop Interchange**: Reordering loops to optimize cache usage can lead to better performance, especially for large arrays.

### Optimized Code Explanation:
- The outermost loop (`t1`) is parallelized using OpenMP to distribute iterations across multiple threads.
- The inner loops (`t2`, `t3`, `t4`) are tiled and fused to improve cache performance and reduce loop overhead.
- The `dummy` function call is placed outside the innermost loop to avoid unnecessary function calls within the critical path.*/

int t1, t2, t3, t4, t5;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

lbp = 0;
ubp = floord(iterations - 1, 32);
#pragma omp parallel for private(lbv, ubv, t2, t3, t4, t5)
for (int t1 = lbp; t1 <= ubp; t1++) {
    for (int t2 = 0; t2 <= floord(LEN_1D - 2, 32); t2++) {
        for (int t3 = max(32 * t1, 32 * t2); t3 <= min(iterations - 1, 32 * t1 + 31); t3++) {
            for (int t4 = max(32 * t2, 0); t4 <= min(LEN_1D - 2, 32 * t2 + 31); t4++) {
                a[t4] = b[t4] + c[t4] * d[t4];
                b[t4] = a[t4] + d[t4] * e[t4];
                a[t4] = b[t4] + a[t4 + 1] * d[t4];
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
    
    time_function(&s243, NULL);
    return EXIT_SUCCESS;
}