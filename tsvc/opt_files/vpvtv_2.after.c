
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


real_t vpvtv(struct args_t * func_args)
{

//    control loops
//    vector plus vector times vector

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
/*### Explanation of the Optimized Code:
1. **Parallelization**: The outer loop (`nl`) is parallelized using OpenMP's `#pragma omp parallel for` directive. This allows multiple threads to execute the iterations of the loop concurrently, which can significantly improve performance on multi-core systems.
2. **Private Variables**: The variables `lbv`, `ubv`, and `i` are declared as private within the OpenMP parallel region to ensure that each thread has its own copy of these variables, avoiding race conditions.
3. **Loop Bounds**: The loop bounds `lbp` and `ubp` are used to define the range of iterations that will be executed in parallel. This is a common technique to manage loop bounds explicitly, which can help in more fine-grained control over the parallel execution.

### Learning from the Examples:
1. **Parallel Loops**: The examples demonstrate the use of OpenMP to parallelize nested loops, which is a common technique for improving performance in multi-threaded environments.
2. **Private Variables**: Declaring loop variables as private within the OpenMP parallel region ensures that each thread operates on its own copy of these variables, preventing data races.
3. **Loop Bounds Management**: Explicitly managing loop bounds (`lbp` and `ubp`) can help in optimizing the distribution of work among threads, especially when dealing with irregular loop structures or dependencies.

By applying these techniques, the optimized code leverages parallelism to enhance the performance of the original loop structure.*/

#include <omp.h>

#define ITERATIONS 100
#define LEN_1D 1000

void dummy(double*, double*, double*, double*, double*, double*, double*, double*, double);

int main() {
    double a[LEN_1D], b[LEN_1D], c[LEN_1D], d[LEN_1D], e[LEN_1D], aa[LEN_1D], bb[LEN_1D], cc[LEN_1D];
    int iterations = ITERATIONS;

    int nl, i;
    int lbp, ubp;
    register int lbv, ubv;

    lbp = 0;
    ubp = iterations - 1;

    #pragma omp parallel for private(lbv, ubv, i)
    for (nl = lbp; nl <= ubp; nl++) {
        for (i = 0; i < LEN_1D; i++) {
            a[i] += b[i] * c[i];
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    return 0;
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
    
    time_function(&vpvtv, NULL);
    return EXIT_SUCCESS;
}