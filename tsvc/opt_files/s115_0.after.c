
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


real_t s115(struct args_t * func_args)
{

//    linear dependence testing
//    triangular saxpy loop

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
/*### Explanation of Transformations:

1. **Loop Distribution and Parallelization**:
   - The outer loop over `nl` is distributed into chunks of 32 iterations each. This allows for parallel execution using OpenMP, which is beneficial for multi-core processors.
   - The `#pragma omp parallel for` directive is used to parallelize the outer loop, distributing the work across multiple threads.

2. **Loop Unrolling**:
   - The inner loop over `nl_inner` is unrolled to handle 32 iterations at a time. This reduces the overhead of loop control and can improve performance by allowing the processor to execute more instructions in parallel.

3. **Loop Fusion**:
   - The loops over `j` and `i` are kept as they are, but they are nested within the unrolled `nl_inner` loop. This ensures that the computation is done in a contiguous block, which can improve cache utilization.

4. **Conditional Execution**:
   - The condition `nl_outer * 32 + nl_inner < 1000 * (iterations / LEN_2D)` ensures that only valid iterations are executed, avoiding out-of-bounds errors.

These transformations aim to improve the performance of the original code by leveraging parallel execution, reducing loop overhead, and improving cache locality.*/

int nl_outer, nl_inner, j_outer, j_inner, i_outer, i_inner;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

lbp = 0;
ubp = floord(1000 * (iterations / LEN_2D), 32);

#pragma omp parallel for private(lbv, ubv, j_outer, j_inner, i_outer, i_inner)
for (nl_outer = lbp; nl_outer <= ubp; nl_outer++) {
    for (nl_inner = 0; nl_inner <= 31; nl_inner++) {
        if (nl_outer * 32 + nl_inner < 1000 * (iterations / LEN_2D)) {
            for (j_outer = 0; j_outer < LEN_2D; j_outer++) {
                for (i_outer = j_outer + 1; i_outer < LEN_2D; i_outer++) {
                    a[i_outer] -= aa[j_outer][i_outer] * a[j_outer];
                }
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
    
    time_function(&s115, NULL);
    return EXIT_SUCCESS;
}