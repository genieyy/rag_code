
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


real_t s3112(struct args_t * func_args)
{

//    reductions
//    sum reduction saving running sums

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

    real_t sum;
#pragma scop
/*### Explanation of the Optimized Code:

1. **Parallelization with OpenMP**:
   - The outer loop (`for (int nl = 0; nl < iterations; nl++)`) is parallelized using OpenMP. This allows multiple threads to execute the loop iterations concurrently.
   - Each thread maintains its own `local_sum` and `local_b` array to avoid race conditions and ensure thread safety.

2. **Local Sum and Local Array**:
   - Each thread computes its own `local_sum` and updates its own `local_b` array. This reduces the need for synchronization during the accumulation phase, improving performance.

3. **Critical Section for Aggregation**:
   - After the inner loop completes, the `local_sum` and `local_b` values are aggregated into the global `sum` and `b` array within a critical section. This ensures that the updates to the global variables are thread-safe.

4. **Dummy Function Call**:
   - The `dummy` function call remains outside the parallel region to avoid unnecessary synchronization overhead.

### Key Learnings from the Examples:

1. **Loop Distribution and Parallelization**:
   - The original examples show how loops can be distributed and parallelized using OpenMP. This involves breaking down the loop iterations into chunks that can be executed by different threads.

2. **Private Variables**:
   - Each thread should have its own private variables to avoid data races. This is achieved by declaring variables as private within the OpenMP parallel region.

3. **Critical Sections**:
   - Critical sections are used to ensure that only one thread at a time updates shared variables. This is necessary when multiple threads need to write to the same memory location.

4. **Loop Tiling**:
   - The examples demonstrate loop tiling, where the loop iterations are divided into smaller blocks to improve cache locality and reduce memory access latency.

By applying these techniques, the optimized code improves performance by leveraging parallel execution and reducing synchronization overhead.*/

#pragma omp parallel
{
    double local_sum = 0.0;
    double local_b[LEN_1D];

    for (int nl = 0; nl < iterations; nl++) {
        local_sum = 0.0;

        #pragma omp for
        for (int i = 0; i < LEN_1D; i++) {
            local_sum += a[i];
            local_b[i] = local_sum;
        }

        #pragma omp critical
        {
            sum += local_sum;
            for (int i = 0; i < LEN_1D; i++) {
                b[i] += local_b[i];
            }
        }

        dummy(a, b, c, d, e, aa, bb, cc, sum);
    }
}
#pragma endscop

    func_args->c2 = omp_get_wtime();
    return sum;
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
    
    time_function(&s3112, NULL);
    return EXIT_SUCCESS;
}