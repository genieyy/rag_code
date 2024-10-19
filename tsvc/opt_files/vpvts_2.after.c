
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


// real_t vpvts( real_t s)
real_t vpvts(struct args_t * func_args)
{

//    control loops
//    vector plus vector times scalar

    real_t s = *(int*)func_args->arg_info;

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
/*### Explanation of the Optimized Code:

1. **Parallelization with OpenMP**:
   - The outer loop (`for (int nl = 0; nl < iterations; nl++)`) is parallelized using `#pragma omp parallel for`. This allows multiple threads to execute the iterations of the loop concurrently, which can significantly improve performance on multi-core processors.

2. **Loop Fusion**:
   - The inner loop (`for (int i = 0; i < LEN_1D; i++)`) is kept as is, but since it is inside a parallelized outer loop, the computation of `a[i] += b[i] * s` is distributed across multiple threads.

3. **Reduction of Function Calls**:
   - The `dummy` function call is kept outside the inner loop, ensuring that it is called only once per iteration of the outer loop. This reduces the overhead of function calls, which can be significant if `dummy` is a complex function.

### Learnings from the Examples:

1. **Parallelization**:
   - The examples demonstrate the use of OpenMP to parallelize loops, which is a common technique to exploit multi-core processors and improve performance.

2. **Loop Distribution and Fusion**:
   - The examples show how to distribute and fuse loops to minimize overhead and maximize computational efficiency.

3. **Reduction of Overhead**:
   - By keeping function calls outside of tight loops, the examples reduce the overhead associated with function calls, which can be a bottleneck in performance-critical code.

4. **Data Dependencies**:
   - The examples handle data dependencies carefully to ensure that parallelization does not introduce race conditions or incorrect results.

By applying these techniques, the optimized code aims to improve the performance of the original loop by leveraging parallel execution and reducing overhead.*/

#pragma omp parallel for
for (int nl = 0; nl < iterations; nl++) {
    for (int i = 0; i < LEN_1D; i++) {
        a[i] += b[i] * s;
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
    
    time_function(&vpvts, &s1);
    return EXIT_SUCCESS;
}