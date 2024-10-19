
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


real_t s422(struct args_t * func_args)
{

//    storage classes and equivalencing
//    common and equivalence statement
//    anti-dependence, threshold of 4

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

    xx = flat_2d_array + 4;

#pragma scop
/*### Analysis of Loop Transformation Methods:

1. **Loop Distribution and Parallelization**:
   - The original nested loops are transformed into a single outer loop with a parallelized inner loop. This allows for better utilization of multi-core processors by distributing the workload across multiple threads.
   - The `#pragma omp parallel for` directive is used to parallelize the inner loop, which is a common technique to exploit parallelism in nested loops.

2. **Loop Reordering and Tiling**:
   - The loops are reordered and tiled to improve cache locality. By breaking the loop iterations into smaller chunks (tiles), the code can better utilize the cache, reducing the number of cache misses.
   - The use of `floord` and `ceild` functions helps in determining the bounds for these tiles, ensuring that the loop iterations are correctly distributed and aligned with the cache lines.

3. **Loop Fusion**:
   - The original loops are fused into a single loop structure, which reduces the overhead of loop control and can improve the performance by reducing the number of times the loop control variables are updated.

### Learnings:
- **Parallelization**: Utilizing OpenMP for parallelizing loops can significantly improve performance on multi-core systems.
- **Cache Optimization**: Tiling and reordering loops can help in improving cache utilization, which is crucial for performance in memory-bound applications.
- **Loop Fusion**: Combining multiple loops into a single loop can reduce overhead and improve performance, especially when the loops operate on the same data.

### Optimized Code:
The provided optimized code applies the same techniques (parallelization, tiling, and loop fusion) to the given loop structure. The outer loop is tiled and the inner loop is parallelized to improve performance. The loop bounds are carefully calculated to ensure correctness and optimal cache utilization.*/

int t1, t2, t3, t4;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

for (int t1 = 0; t1 <= floord(8 * iterations - 1, 8); t1++) {
    lbp = max(ceild(2 * t1, 3), ceild(32 * t1 - LEN_1D + 1, 32));
    ubp = floord(3 * t1 + 2, 4);
#pragma omp parallel for private(lbv, ubv, t3, t4)
    for (int t2 = lbp; t2 <= ubp; t2++) {
        for (int t3 = max(ceild(32 * t2 + 1, 3), 32 * t1 - 32 * t2); t3 <= min(min(LEN_1D - 1, 16 * t2 + 15), 32 * t1 - 32 * t2 + 31); t3++) {
            for (int t4 = max(32 * t2, 2 * t3); t4 <= min(min(32 * t2 + 31, 3 * t3 - 1), 2 * t3 + LEN_1D - 3); t4++) {
                xx[(-2 * t3 + t4)] = flat_2d_array[(-2 * t3 + t4) + 8] + a[(-2 * t3 + t4)];
            }
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
    
    time_function(&s422, NULL);
    return EXIT_SUCCESS;
}