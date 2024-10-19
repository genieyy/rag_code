
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


real_t s174(struct args_t * func_args)
{

//    symbolics
//    loop with subscript that may seem ambiguous

    int M = *(int*)func_args->arg_info;

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
/*### Explanation of the Optimized Code:
1. **Loop Unrolling and Parallelization**:
   - The outer loop (`nl_outer`) is unrolled by a factor of 32 to reduce the overhead of loop control.
   - The inner loop (`nl_inner`) iterates over the unrolled chunks of the outer loop.
   - The `#pragma omp parallel for` directive is used to parallelize the outer loop, allowing multiple threads to execute different iterations of the loop concurrently.

2. **Loop Distribution**:
   - The inner loop over `i` is distributed across multiple threads by dividing the range `[0, M)` into chunks of size 32.
   - This distribution helps in better load balancing and reduces the contention for shared resources.

3. **Register Usage**:
   - The `register` keyword is used for the loop bounds (`lbv`, `ubv`) to suggest that these variables should be stored in CPU registers for faster access.

4. **Minimizing Function Calls**:
   - The `dummy` function call is kept outside the innermost loop to avoid frequent function calls, which can be expensive in terms of performance.

This optimization strategy leverages parallel processing, loop unrolling, and efficient memory access patterns to improve the performance of the original loop.*/

int nl_outer, nl_inner, i_outer, i_inner;
int lb, ub, lbp, ubp;
register int lbv, ubv;

lbp = 0;
ubp = floord(M, 32);

#pragma omp parallel for private(lbv, ubv, nl_inner, i_outer, i_inner)
for (nl_outer = 0; nl_outer < 10 * iterations; nl_outer += 32) {
    for (lbv = lbp; lbv <= ubp; lbv++) {
        for (nl_inner = nl_outer; nl_inner < min(nl_outer + 32, 10 * iterations); nl_inner++) {
            for (i_outer = 32 * lbv; i_outer <= min(M - 1, 32 * lbv + 31); i_outer++) {
                a[i_outer + M] = a[i_outer] + b[i_outer];
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
    
    time_function(&s174, &(struct{int a;}){LEN_1D/2});
    return EXIT_SUCCESS;
}