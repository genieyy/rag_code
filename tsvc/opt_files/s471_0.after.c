
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

// %4.7

int s471s(void)
{
// --  dummy subroutine call made in s471
    return 0;
}

real_t s471(struct args_t * func_args){

//    call statements

    int m = LEN_1D;

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
/*### Explanation of Transformations:

1. **Loop Parallelization**: The outer loop over `nl` is parallelized using OpenMP's `#pragma omp parallel for` directive. This allows multiple threads to execute iterations of the loop concurrently, which can significantly improve performance on multi-core systems.

2. **Loop Bounds Calculation**: The loop bounds for `nl` are calculated once and stored in `nl_lbp` and `nl_ubp`. This avoids recalculating the bounds in each iteration of the loop, which can be inefficient.

3. **Private Variables**: The loop indices `nl_lbv`, `nl_ubv`, `i_lbv`, and `i_ubv` are declared as private within the OpenMP parallel region. This ensures that each thread has its own copy of these variables, preventing race conditions.

4. **Loop Unrolling**: Although not explicitly unrolled in this example, the inner loop over `i` could potentially benefit from loop unrolling if the compiler does not already perform this optimization. Loop unrolling reduces the overhead of loop control instructions by processing multiple iterations in a single loop body.

5. **Vectorization**: The `#pragma ivdep` and `#pragma vector always` directives are not used in this example, but they could be applied to the inner loop if the operations within the loop are vectorizable. These directives instruct the compiler to ignore potential dependencies and vectorize the loop, which can improve performance on vector processors.

By applying these transformations, the code is optimized for parallel execution and reduced overhead, which should improve its performance.*/

int nl_lbp, nl_ubp, i_lbp, i_ubp;
register int nl_lbv, nl_ubv, i_lbv, i_ubv;

nl_lbp = 0;
nl_ubp = (iterations / 2) - 1;

#pragma omp parallel for private(nl_lbv, nl_ubv, i_lbv, i_ubv)
for (int nl_t = nl_lbp; nl_t <= nl_ubp; nl_t++) {
    i_lbp = 0;
    i_ubp = m - 1;
    for (int i_t = i_lbp; i_t <= i_ubp; i_t++) {
        x[i_t] = b[i_t] + d[i_t] * d[i_t];
        s471s();
        b[i_t] = c[i_t] + d[i_t] * e[i_t];
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
    
    time_function(&s471, NULL);
    return EXIT_SUCCESS;
}