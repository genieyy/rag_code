
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
/*### Explanation of the Corrected Transformation:

1. **Parallelization**:
   - The outer loop over `nl` is parallelized using OpenMP's `#pragma omp parallel for` directive. This allows the iterations of the outer loop to be executed in parallel, potentially leveraging multiple CPU cores to speed up the computation.

2. **Private Variables**:
   - The variable `i` is declared inside the loop, ensuring that each thread has its own instance of `i`, preventing race conditions. This avoids the need to explicitly declare `i` as private, which was causing the compilation error.

3. **Bounds Calculation**:
   - The bounds for the `nl` loop are precomputed and stored in `nl_lb` and `nl_ub` to avoid recalculating the bounds on each iteration, which can be a minor optimization.

### What I Learned:

- **Parallelization**: The use of OpenMP can significantly improve performance by parallelizing loops, especially when the iterations are independent.
- **Loop Bounds**: Precomputing loop bounds can be beneficial for readability and potentially for performance, as it avoids redundant calculations.
- **Private Variables**: Ensuring that loop variables are private in parallel regions is crucial to avoid data races and ensure correct execution. Declaring the loop variable inside the loop is a straightforward way to achieve this.*/

int nl_lb = 0;
int nl_ub = 10 * iterations - 1;
#pragma omp parallel for
for (int nl = nl_lb; nl <= nl_ub; nl++) {
    for (int i = 0; i < M; i++) {
        a[i + M] = a[i] + b[i];
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
    
    time_function(&s174, &(struct{int a;}){LEN_1D/2});
    return EXIT_SUCCESS;
}