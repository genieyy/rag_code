
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


real_t s317(struct args_t * func_args)
{

//    reductions
//    product reductio vectorize with
//    1. scalar expansion of factor, and product reduction
//    2. closed form solution: q = factor**n

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

    real_t q;
#pragma scop
/*### Explanation of the Optimized Code:

1. **Parallelization**: The outer loop is parallelized using OpenMP (`#pragma omp parallel for`). This allows multiple threads to execute the loop iterations concurrently, which can significantly improve performance on multi-core processors.

2. **Private Variable**: The variable `q` is declared as private within the parallel region. This ensures that each thread has its own separate instance of `q`, preventing race conditions and ensuring thread safety.

3. **Loop Unrolling**: The inner loop is not unrolled in this example, but if `LEN_1D` is known to be a multiple of a small number (e.g., 4 or 8), loop unrolling could be considered to reduce the number of iterations and improve performance. However, this optimization is not applied here to keep the code simple and maintainable.

4. **Vectorization**: The inner loop could potentially benefit from vectorization, but this is typically handled by the compiler if the loop is simple and the hardware supports it. The `#pragma ivdep` and `#pragma vector always` directives are not used here, as they are more relevant for more complex loops with dependencies.

By applying these optimizations, the code is expected to run faster, especially on systems with multiple CPU cores.*/

#pragma omp parallel for private(q)
for (int nl = 0; nl < 5 * iterations; nl++) {
    double q = 1.0;
    for (int i = 0; i < LEN_1D / 2; i++) {
        q *= 0.99;
    }
    dummy(a, b, c, d, e, aa, bb, cc, q);
}
#pragma endscop

    func_args->c2 = omp_get_wtime();
    return q;
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
    
    time_function(&s317, NULL);
    return EXIT_SUCCESS;
}