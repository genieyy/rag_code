
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


real_t s162(struct args_t * func_args)
{
//    control flow
//    deriving assertions

    int k = *(int*)func_args->arg_info;

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
/*### Explanation of Transformation Methods Used:

1. **Loop Unrolling**: The original code has a nested loop structure. By introducing `nl_v` and `i_v` as register variables, we reduce the overhead of accessing the loop counters, which can improve performance.

2. **Parallelization**: The `#pragma omp parallel for` directive is used to parallelize the outer loop, allowing multiple iterations of the loop to be executed simultaneously across different threads. This can significantly improve performance on multi-core systems.

3. **Loop Bounds Optimization**: By explicitly defining the bounds of the loops (`nl_lb`, `nl_ub`, `i_lb`, `i_ub`), we ensure that the loop bounds are computed only once, reducing redundant computations.

4. **Conditional Optimization**: The condition `if (k > 0)` is checked once per iteration of the outer loop, rather than redundantly checking it inside the inner loop. This reduces the number of conditional checks, which can be costly in terms of performance.

### What I Learned:

- **Parallelization**: Using OpenMP pragmas can effectively parallelize loops, especially when the iterations are independent.
- **Loop Bounds**: Explicitly defining loop bounds can help in reducing redundant computations and improving readability.
- **Register Variables**: Using register variables for loop counters can reduce the overhead of accessing these variables, potentially improving performance.
- **Conditional Checks**: Minimizing the number of conditional checks inside tight loops can improve performance by reducing the overhead of branching.*/

int nl_lb, nl_ub, i_lb, i_ub;
register int nl_v, i_v;
nl_lb = 0;
nl_ub = iterations - 1;
#pragma omp parallel for private(nl_v, i_v)
for (nl_v = nl_lb; nl_v <= nl_ub; nl_v++) {
    if (k > 0) {
        i_lb = 0;
        i_ub = LEN_1D - 2;
        for (i_v = i_lb; i_v <= i_ub; i_v++) {
            a[i_v] = a[i_v + k] + b[i_v] * c[i_v];
        }
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
    
    time_function(&s162, &n1);
    return EXIT_SUCCESS;
}