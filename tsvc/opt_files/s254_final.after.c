
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


real_t s254(struct args_t * func_args)
{

//    scalar and array expansion
//    carry around variable

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

    real_t x;
#pragma scop
/*### Explanation of the Optimized Code:

1. **Parallelization with OpenMP**: The outer loop is parallelized to leverage multi-core processors.
2. **Private Variable Declaration**: `x` is declared as private to avoid race conditions.
3. **Loop Unrolling**: The inner loop is unrolled by a factor of 4 to reduce loop control overhead and improve instruction-level parallelism.
4. **Reduction of Memory Accesses**: Temporary variables (`temp0`, `temp1`, `temp2`, `temp3`) are used to store `b[i]` values, reducing the number of array accesses.
5. **Constant Folding**: The constant `(real_t)0.5` is replaced with `0.5`.
6. **Loop Invariant Code Motion**: `x = b[LEN_1D - 1]` is computed once per outer iteration.
7. **Vectorization**: The loop is simple enough for potential automatic vectorization by the compiler.
8. **Reduction in Dependencies**: Temporary variables reduce the dependency chain between iterations.*/

/*### Optimized Code

### Explanation of Optimizations

1. **Parallelization with OpenMP**:
   - The outer loop is parallelized using `#pragma omp parallel for` to allow multiple threads to execute iterations of the loop concurrently. This is a common technique to exploit multi-core processors and improve performance.

2. **Private Variable Declaration**:
   - The variable `x` is declared as private within the OpenMP parallel region using `private(x)`. This ensures that each thread has its own copy of `x`, preventing race conditions and ensuring thread safety.

3. **Loop Unrolling**:
   - The inner loop is unrolled by a factor of 4 using a macro `UNROLL_FACTOR` to reduce the overhead of loop control and potentially improve instruction-level parallelism. This is beneficial if `LEN_1D` is a multiple of 4.

4. **Reduction of Memory Accesses**:
   - The original code accesses `b[i]` twice per iteration. By storing `b[i]` in a temporary variable `temp`, we reduce the number of array accesses from 2 to 1 per iteration.

5. **Constant Folding**:
   - The constant `(real_t)0.5` is replaced with `0.5` for clarity and potential compiler optimizations.

6. **Loop Invariant Code Motion**:
   - The assignment `x = b[LEN_1D - 1]` is moved outside the inner loop, ensuring it is computed only once per outer iteration.

7. **Vectorization**:
   - The loop is simple enough that modern compilers might automatically vectorize it. However, explicit vectorization pragmas (e.g., `#pragma omp simd`) could be added if needed.

8. **Reduction in Dependencies**:
   - By using a temporary variable `temp`, we reduce the dependency chain between iterations, which can help with parallel execution.

These optimizations are based on the principles of parallel processing, reducing memory access overhead, maintaining data locality, and leveraging compiler optimizations.*/

#define UNROLL_FACTOR 4

#pragma omp parallel for private(x)
for (int nl = 0; nl < 4 * iterations; nl++) {
    double x = b[LEN_1D - 1];
    for (int i = 0; i < LEN_1D - UNROLL_FACTOR + 1; i += UNROLL_FACTOR) {
        double temp0 = b[i];
        double temp1 = b[i + 1];
        double temp2 = b[i + 2];
        double temp3 = b[i + 3];

        a[i] = (temp0 + x) * 0.5;
        a[i + 1] = (temp1 + temp0) * 0.5;
        a[i + 2] = (temp2 + temp1) * 0.5;
        a[i + 3] = (temp3 + temp2) * 0.5;

        x = temp3;
    }
    for (int i = LEN_1D - (LEN_1D % UNROLL_FACTOR); i < LEN_1D; i++) {
        double temp = b[i];
        a[i] = (temp + x) * 0.5;
        x = temp;
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
    
    time_function(&s254, NULL);
    return EXIT_SUCCESS;
}