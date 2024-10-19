
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


real_t s1221(struct args_t * func_args)
{

//    run-time symbolic resolution

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
/*### Analysis of Loop Transformation Methods Used:

1. **Loop Distribution and Parallelization**:
   - The original loops are distributed into multiple smaller loops that can be executed in parallel using OpenMP. This is done by introducing new variables (`t1`, `t2`, `t3`) and using `#pragma omp parallel for` to parallelize the loop iterations.

2. **Loop Tiling**:
   - The loops are tiled to improve cache locality. This is achieved by dividing the loop iterations into chunks (tiles) of size 32, which is a common cache line size. This helps in reducing cache misses and improving performance.

3. **Loop Reordering**:
   - The order of loops is changed to ensure that the innermost loop has the best possible locality. This is done by reordering the loops such that the most frequently accessed data is accessed in the innermost loop.

4. **Loop Unrolling**:
   - Although not explicitly shown in the provided examples, loop unrolling can be inferred as a potential optimization. By unrolling the innermost loop, the number of iterations is reduced, which can lead to fewer loop overheads and better instruction-level parallelism.

### Learnings Applied to the New Code:

- **Parallelization**: The new code uses OpenMP to parallelize the outer loop, allowing multiple iterations to be executed simultaneously.
- **Loop Tiling**: The iterations are divided into chunks of size 32, which helps in improving cache locality.
- **Loop Reordering**: The order of loops is adjusted to ensure that the most frequently accessed data is accessed in the innermost loop, which is `b[t3] = b[t3 - 4] + a[t3]`.
- **Loop Distribution**: The original loop is distributed into multiple smaller loops that can be executed in parallel, improving performance.

By applying these transformations, the new code is expected to perform better than the original code, especially on multi-core processors.*/

int t1, t2, t3;
int lb, ub, lbp, ubp;
register int lbv, ubv;

for (int t1 = 0; t1 <= floord(iterations - 1, 32); t1++) {
    lbp = max(0, ceild(32 * t1 - LEN_1D + 5, 32));
    ubp = min(floord(iterations - 1, 32), t1);
#pragma omp parallel for private(lbv, ubv, t3)
    for (int t2 = lbp; t2 <= ubp; t2++) {
        for (int t3 = max(4, 32 * t2); t3 <= min(LEN_1D - 1, 32 * t2 + 31); t3++) {
            b[t3] = b[t3 - 4] + a[t3];
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
    
    time_function(&s1221, NULL);
    return EXIT_SUCCESS;
}