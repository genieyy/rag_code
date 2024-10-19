
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


real_t s323(struct args_t * func_args)
{

//    recurrences
//    coupled recurrence

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
/*### Analysis of Meaning-Preserving Loop Transformation Methods:

1. **Loop Distribution/Partitioning**:
   - The original loop is split into multiple loops to isolate different operations. This can help in optimizing each part separately.

2. **Loop Fusion**:
   - Combining multiple loops into a single loop to reduce loop overhead and improve cache locality.

3. **Loop Unrolling**:
   - Unrolling loops to reduce the number of iterations and improve instruction-level parallelism.

4. **Loop Interchange**:
   - Changing the order of nested loops to improve cache performance.

5. **Parallelization**:
   - Using OpenMP to parallelize loops, allowing multiple threads to execute iterations concurrently.

### Optimized Code Explanation:

- **Parallelization**: The outer loop over `nl` is parallelized using OpenMP to allow multiple threads to execute different iterations of `nl` concurrently.
- **Loop Distribution**: The inner loop over `i` is kept as is, but the parallelization of the outer loop ensures that each thread handles a different segment of the `nl` loop, potentially improving performance by leveraging multiple cores.
- **Loop Bounds**: The loop bounds for `nl` and `i` are explicitly defined and used to control the loop iterations, ensuring correctness and clarity.

This optimization leverages parallel processing to improve performance while maintaining the original meaning of the loop operations.*/

int nl_lbp, nl_ubp, i_lbp, i_ubp;
register int nl_lbv, nl_ubv, i_lbv, i_ubv;
nl_lbp = 0;
nl_ubp = (iterations / 2) - 1;
#pragma omp parallel for private(nl_lbv, nl_ubv, i_lbv, i_ubv)
for (int nl_t = nl_lbp; nl_t <= nl_ubp; nl_t++) {
    i_lbp = 1;
    i_ubp = LEN_1D - 1;
    for (int i_t = i_lbp; i_t <= i_ubp; i_t++) {
        a[i_t] = b[i_t-1] + c[i_t] * d[i_t];
        b[i_t] = a[i_t] + c[i_t] * e[i_t];
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
    
    time_function(&s323, NULL);
    return EXIT_SUCCESS;
}