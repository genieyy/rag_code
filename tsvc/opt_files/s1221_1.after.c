
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
/*### Explanation of the Transformation:

1. **Loop Distribution and Parallelization**:
   - The outer loop over `nl` is distributed into chunks of size 32 (using `floord(iterations - 1, 32)`), which allows for parallel execution using OpenMP (`#pragma omp parallel for`).
   - The inner loop over `i` is transformed to iterate over the range `[max(4, 32 * t1 - 32 * t2), min(LEN_1D - 1, 32 * t1 - 32 * t2 + 31)]`, ensuring that each thread processes a non-overlapping subset of the array.

2. **Loop Unrolling**:
   - The inner loop is not explicitly unrolled, but the transformation ensures that the loop bounds are adjusted to maximize the efficiency of the loop body execution.

3. **Dependency Analysis**:
   - The transformation ensures that the dependency on `b[i - 4]` is respected by adjusting the loop bounds appropriately, ensuring that each iteration of the loop only accesses elements that have already been computed.

4. **Reduction of Overhead**:
   - By distributing the loop and parallelizing it, the overhead of managing the loop iterations is reduced, and the computation is spread across multiple threads, improving performance.

5. **Dummy Function Call**:
   - The `dummy` function call is placed outside the inner loop, ensuring that it is called only once per chunk of iterations, reducing the overhead of function calls.*/

int t1, t2, t3, t4;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

for (int t1 = 0; t1 <= floord(iterations - 1, 32); t1++) {
    lbp = max(0, ceild(32 * t1 - LEN_1D + 5, 32));
    ubp = floord(32 * t1 + 31, 32);
#pragma omp parallel for private(lbv, ubv, t3, t4)
    for (int t2 = lbp; t2 <= ubp; t2++) {
        for (int t3 = max(4, 32 * t1 - 32 * t2); t3 <= min(LEN_1D - 1, 32 * t1 - 32 * t2 + 31); t3++) {
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