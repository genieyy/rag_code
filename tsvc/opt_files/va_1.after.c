
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


real_t va(struct args_t * func_args)
{

//    control loops
//    vector assignment

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
/*### Explanation of Loop Transformation Methods Used

1. **Loop Unrolling**: The outer loop is unrolled by a factor of 32 (similar to the examples provided). This reduces the number of iterations of the outer loop, which can be beneficial for performance, especially when combined with parallelization.

2. **Parallelization**: The outer loop is parallelized using OpenMP (`#pragma omp parallel for`). This allows multiple threads to execute the loop in parallel, which can significantly improve performance on multi-core systems.

3. **Loop Distribution**: The inner loop is distributed across the unrolled iterations of the outer loop. This ensures that each iteration of the outer loop handles a chunk of the inner loop's work, which can be more efficient than a single large loop.

4. **Conditional Execution**: The condition `if (nl_outer * 32 + nl_inner < iterations * 10)` ensures that only valid iterations are executed, preventing out-of-bounds access and ensuring correctness.

5. **Variable Scope**: The variable `i` is declared inside the inner loop to ensure it is in scope for the loop body, avoiding the compilation error.

These transformations are inspired by the examples provided, where similar techniques were used to optimize nested loops. The goal is to reduce the overhead of loop control and maximize the potential for parallel execution.*/

int nl_outer, nl_inner;
int lb, ub, lbp, ubp;
register int lbv, ubv;

lbp = 0;
ubp = floord(iterations * 10, 32);
#pragma omp parallel for private(lbv, ubv, nl_inner, i)
for (nl_outer = lbp; nl_outer <= ubp; nl_outer++) {
    for (nl_inner = 0; nl_inner <= 31; nl_inner++) {
        if (nl_outer * 32 + nl_inner < iterations * 10) {
            int i; // Declare i here to ensure it is in scope
            for (i = 0; i < LEN_1D; i++) {
                a[i] = b[i];
            }
            dummy(a, b, c, d, e, aa, bb, cc, 0.);
        }
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
    
    time_function(&va, NULL);
    return EXIT_SUCCESS;
}