
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


real_t s313(struct args_t * func_args)
{

//    reductions
//    dot product

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

    real_t dot;
#pragma scop
/*### Optimized Code



### Analysis of Loop Transformation Methods

1. **Loop Unrolling**: The original code can be unrolled to reduce the overhead of loop control. However, in this case, the loop is already simple and unrolling might not provide significant benefits.

2. **Vectorization**: The inner loop can be vectorized to take advantage of SIMD instructions. This is done using compiler directives like `#pragma ivdep` and `#pragma vector always` in the provided examples.

3. **Parallelization**: The outer loop can be parallelized using OpenMP to distribute the workload across multiple threads. This is done using `#pragma omp parallel for` in the provided examples.

4. **Reduction**: The inner loop accumulates a sum, which can be optimized using a reduction clause in OpenMP. However, since the outer loop is not parallelized, this is not applicable here.

5. **Loop Fusion**: If possible, combining multiple loops into one can reduce overhead. However, in this case, the loops are independent and cannot be fused.

### Application of Methods

- **Vectorization**: The inner loop is vectorized using a temporary variable `dot_sum` to accumulate the sum. This reduces the overhead of repeatedly accessing and updating the `dot` variable.
- **Parallelization**: The outer loop is not parallelized here, but it could be parallelized using OpenMP if `dummy` is thread-safe.

### Optimized Code Explanation

- **Temporary Variable**: A temporary variable `dot_sum` is used to accumulate the sum in the inner loop. This reduces the overhead of repeatedly accessing and updating the `dot` variable.
- **Initialization**: The `dot_sum` variable is initialized to `0.0` at the start of each iteration of the outer loop to ensure correct accumulation.
- **Function Call**: The `dummy` function is called with the accumulated `dot_sum` value.

This optimization focuses on reducing the overhead of accumulation and ensuring that the loop remains simple and efficient.*/

double dot_sum = 0.0;
for (int nl = 0; nl < iterations * 5; nl++) {
    dot_sum = 0.0;
    for (int i = 0; i < LEN_1D; i++) {
        dot_sum += a[i] * b[i];
    }
    dummy(a, b, c, d, e, aa, bb, cc, dot_sum);
}
#pragma endscop

    func_args->c2 = omp_get_wtime();
    return dot;
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
    
    time_function(&s313, NULL);
    return EXIT_SUCCESS;
}