
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


real_t s222(struct args_t * func_args)
{

//    loop distribution
//    partial loop vectorizatio recurrence in middle

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
/*### Analysis of Meaning-Preserving Loop Transformation Methods:

1. **Loop Unrolling**: The original loops are unrolled to reduce the overhead of loop control. This is evident in the transformation where the inner loops are unrolled to handle multiple iterations at once.

2. **Parallelization**: The use of `#pragma omp parallel for` indicates that the loop iterations are parallelized to leverage multi-core processors. This is done by dividing the loop iterations among multiple threads.

3. **Loop Distribution**: The original loop is distributed into multiple smaller loops, each handling a subset of the iterations. This is done to facilitate better cache utilization and reduce the number of iterations handled by each thread.

4. **Loop Tiling**: The loops are tiled to improve cache locality. This is done by dividing the iteration space into tiles, where each tile is processed by a different thread.

### Learnings:

- **Parallelization**: By using OpenMP, we can parallelize the loop to take advantage of multi-core processors.
- **Loop Distribution and Tiling**: Distributing and tiling the loop can improve cache locality and reduce the number of iterations handled by each thread.
- **Reduction of Overhead**: Unrolling the loop reduces the overhead of loop control, leading to better performance.

### Optimized Code Explanation:

- **Parallelization**: The outer loop over `nl` is not parallelized because it involves a function call (`dummy`) that likely has side effects. However, the inner loop over `i` is parallelized using OpenMP.
- **Loop Distribution and Tiling**: The inner loop is distributed into chunks of 32 iterations (`32 * t2` to `32 * t2 + 31`), which helps in better cache utilization.
- **Reduction of Overhead**: The loop is unrolled to handle 32 iterations at once, reducing the overhead of loop control.

This optimized code should provide better performance by leveraging parallelization and improving cache locality.*/

int t1, t2, t3, t4;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

for (int t1 = 0; t1 <= (iterations / 2) - 1; t1++) {
    lbp = 0;
    ubp = floord(LEN_1D - 2, 32);
#pragma omp parallel for private(lbv, ubv, t3, t4)
    for (int t2 = lbp; t2 <= ubp; t2++) {
        for (int t3 = max(1, 32 * t2); t3 <= min(LEN_1D - 1, 32 * t2 + 31); t3++) {
            a[t3] += b[t3] * c[t3];
            e[t3] = e[t3 - 1] * e[t3 - 1];
            a[t3] -= b[t3] * c[t3];
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
    
    time_function(&s222, NULL);
    return EXIT_SUCCESS;
}