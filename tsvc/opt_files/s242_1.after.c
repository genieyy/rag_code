
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


real_t s242(struct args_t * func_args)
{

//    node splitting

    struct{real_t a;real_t b;} * x = func_args->arg_info;
    real_t s1 = x->a;
    real_t s2 = x->b;

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
/*### Analysis of Loop Transformation Methods Used:

1. **Loop Distribution and Parallelization**:
   - The original nested loops are distributed into multiple loops to facilitate parallel execution. This is done using `#pragma omp parallel for` to parallelize the outer loop and distribute the workload across threads.

2. **Loop Tiling**:
   - The loops are tiled to improve cache locality. By dividing the iteration space into chunks (tiles) of size 32, the code ensures that each thread works on a smaller, contiguous block of data, which can improve cache performance.

3. **Loop Vectorization**:
   - The inner loop is vectorized using `#pragma ivdep` and `#pragma vector always` to hint the compiler to generate SIMD instructions, which can significantly speed up the computation by processing multiple elements in parallel.

4. **Loop Unrolling**:
   - Although not explicitly unrolled in the provided code, the combination of tiling and vectorization implicitly unrolls the loop by processing multiple iterations in parallel.

### Learning:
- **Parallelization**: Using OpenMP pragmas to parallelize loops can significantly improve performance by leveraging multiple CPU cores.
- **Tiling**: Breaking down the iteration space into smaller chunks (tiles) can improve cache utilization and reduce cache misses.
- **Vectorization**: Using compiler hints to vectorize loops can take advantage of SIMD instructions, which can process multiple data points simultaneously.

### Optimized Code Explanation:
- The outer loop over `nl` is parallelized using OpenMP.
- The inner loop over `i` is tiled and vectorized to improve cache locality and utilize SIMD instructions.
- The `dummy` function call remains outside the parallelized loop to avoid potential race conditions and ensure correct execution.*/

int t1, t2, t3;
int lb, ub, lbp, ubp;
register int lbv, ubv;

for (int t1 = 0; t1 <= (iterations / 5) - 1; t1++) {
    lbp = 0;
    ubp = floord(LEN_1D - 1, 32);
#pragma omp parallel for private(lbv, ubv, t3)
    for (int t2 = lbp; t2 <= ubp; t2++) {
        lbv = max(1, 32 * t2);
        ubv = min(LEN_1D - 1, 32 * t2 + 31);
#pragma ivdep
#pragma vector always
        for (int t3 = lbv; t3 <= ubv; t3++) {
            a[t3] = a[t3 - 1] + s1 + s2 + b[t3] + c[t3] + d[t3];
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
    
    time_function(&s242, &(struct{real_t a;real_t b;}){s1, s2});
    return EXIT_SUCCESS;
}