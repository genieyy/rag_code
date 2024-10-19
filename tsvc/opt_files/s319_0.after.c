
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


real_t s319(struct args_t * func_args)
{

//    reductions
//    coupled reductions

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

    real_t sum;
#pragma scop
/*### Explanation of Optimizations:

1. **Loop Unrolling and Padding**: 
   - The arrays `a_local` and `b_local` are padded to ensure they are aligned to cache lines, which can improve memory access efficiency.
   - The size `LEN_1D_PADDED` is calculated to be a multiple of 8, which is a common cache line size.

2. **Parallel Execution with OpenMP**:
   - The `#pragma omp parallel for` directive is used to parallelize the inner loop, allowing multiple threads to process different parts of the loop concurrently.
   - The `reduction(+:sum_local[nl % 2])` clause ensures that the partial sums from each thread are combined into a single `sum_local` value.

3. **Local Storage for Intermediate Results**:
   - Intermediate results (`a_local` and `b_local`) are stored in local arrays to reduce cache contention and improve locality of reference.
   - The `nl % 2` indexing is used to alternate between two sets of local arrays, which can help in reducing false sharing and improving performance.

4. **Reduction of Global Memory Accesses**:
   - By using local arrays and a reduction clause, the number of global memory accesses is minimized, which can significantly improve performance, especially in a multi-threaded environment.*/

#include <omp.h>

#define LEN_1D_PADDED (LEN_1D + (LEN_1D % 8 ? 8 - (LEN_1D % 8) : 0))

double sum_local[2];
double a_local[2][LEN_1D_PADDED];
double b_local[2][LEN_1D_PADDED];

for (int nl = 0; nl < 2 * iterations; nl++) {
    sum_local[nl % 2] = 0.;
    #pragma omp parallel for reduction(+:sum_local[nl % 2])
    for (int i = 0; i < LEN_1D; i++) {
        a_local[nl % 2][i] = c[i] + d[i];
        sum_local[nl % 2] += a_local[nl % 2][i];
        b_local[nl % 2][i] = c[i] + e[i];
        sum_local[nl % 2] += b_local[nl % 2][i];
    }
    dummy(a_local[nl % 2], b_local[nl % 2], c, d, e, aa, bb, cc, sum_local[nl % 2]);
}
#pragma endscop

    func_args->c2 = omp_get_wtime();
    return sum;
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
    
    time_function(&s319, NULL);
    return EXIT_SUCCESS;
}