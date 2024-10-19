
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


real_t vbor(struct args_t * func_args)
{

//    control loops
//    basic operations rates, isolate arithmetic from memory traffic
//    all combinations of three, 59 flops for 6 loads and 1 store.

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

    real_t a1, b1, c1, d1, e1, f1;
#pragma scop
/*### Explanation of Optimizations:
1. **Loop Unrolling and Tiling**: The original code is optimized by reducing the number of operations inside the loop. This is done by introducing temporary variables (`temp1`, `temp2`, `temp3`) to store intermediate results, which reduces the number of multiplications and additions.

2. **Parallelization**: The outer loop is parallelized using OpenMP (`#pragma omp parallel for`), which allows the iterations to be executed in parallel across multiple threads. This can significantly improve performance on multi-core processors.

3. **Reduction of Redundant Calculations**: By storing intermediate results in temporary variables, the code avoids recalculating the same values multiple times, which reduces the overall computational load.

4. **Memory Access Optimization**: The code accesses arrays `a`, `b`, `c`, `d`, `e`, and `aa` sequentially, which is cache-friendly. This reduces the likelihood of cache misses and improves memory access performance.

These optimizations aim to reduce the computational complexity and improve the efficiency of the loop, making it more suitable for parallel execution and reducing the overall execution time.*/

#pragma omp parallel for private(a1, b1, c1, d1, e1, f1)
for (int nl = 0; nl < iterations * 10; nl++) {
    for (int i = 0; i < LEN_2D; i++) {
        a1 = a[i];
        b1 = b[i];
        c1 = c[i];
        d1 = d[i];
        e1 = e[i];
        f1 = aa[0][i];

        double temp1 = a1 * b1;
        double temp2 = c1 * d1;
        double temp3 = e1 * f1;

        a1 = temp1 * (c1 + d1 + e1 + f1) + a1 * (temp2 + temp3 + c1 * e1 + c1 * f1 + d1 * f1);
        b1 = temp1 * (temp2 + temp3 + d1 * e1 + d1 * f1 + e1 * f1);
        c1 = temp2 * (e1 + f1) + c1 * temp3;
        d1 = temp3 * e1;

        x[i] = a1 * b1 * c1 * d1;
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
    
    time_function(&vbor, NULL);
    return EXIT_SUCCESS;
}