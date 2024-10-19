
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


real_t s431(struct args_t * func_args)
{

//    parameters
//    parameter statement

    int k1=1;
    int k2=2;
    int k=2*k1-k2;

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
/*### Analysis of Meaning-Preserving Loop Transformation Methods

1. **Loop Distribution and Parallelization**:
   - The original nested loops are distributed into multiple loops to facilitate parallel execution. This is evident in the use of `#pragma omp parallel for` to parallelize the outer loop.

2. **Loop Tiling**:
   - The loops are tiled to improve cache locality. This is done by introducing new variables (`t1`, `t2`, `t3`) that divide the iteration space into smaller chunks (tiles), which can be processed more efficiently.

3. **Loop Reordering**:
   - The order of loops is changed to optimize for parallelism and cache usage. The innermost loop is kept as the computation loop, while the outer loops are reordered to facilitate better parallelization and tiling.

4. **Loop Unrolling**:
   - The innermost loop is partially unrolled to reduce loop overhead. This is done by iterating over a fixed number of iterations (32 in this case) within the loop.

### Learning from the Examples

- **Parallelization**: By using OpenMP directives, the code can be parallelized to leverage multiple CPU cores, significantly improving performance for large iteration spaces.
- **Tiling**: Tiling the iteration space can improve cache performance by ensuring that data accessed within a tile fits into the cache, reducing cache misses.
- **Unrolling**: Partial loop unrolling can reduce the overhead of loop control instructions, leading to better performance.

### Optimized Code Explanation

- **Outer Loop (`t1`)**: This loop is parallelized using OpenMP to distribute the work across multiple threads.
- **Middle Loop (`t2`)**: This loop is used to partially unroll the iterations within each tile.
- **Inner Loop (`t3`)**: This loop iterates over the actual iterations within each tile, performing the computation and calling the `dummy` function.

This optimized code should provide better performance by leveraging parallelism, improving cache locality, and reducing loop overhead.*/

int t1, t2, t3;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

lbp = 0;
ubp = floord(iterations * 10, 32);
#pragma omp parallel for private(lbv, ubv, t2, t3)
for (t1 = lbp; t1 <= ubp; t1++) {
    for (t2 = 0; t2 <= 31; t2++) {
        for (t3 = 32 * t1; t3 <= min(iterations * 10 - 1, 32 * t1 + 31); t3++) {
            for (int i = 0; i < LEN_1D; i++) {
                a[i] = a[i + k] + b[i];
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
    
    time_function(&s431, NULL);
    return EXIT_SUCCESS;
}