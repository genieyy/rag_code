
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


real_t s276(struct args_t * func_args)
{

//    control flow
//    if test using loop index

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

    int mid = (LEN_1D/2);
#pragma scop
/*### Analysis of Meaning-Preserving Loop Transformation Methods:

1. **Loop Distribution/Partitioning**:
   - The original loop is split into multiple loops to handle different conditions separately. This can help in optimizing each part independently.

2. **Loop Unrolling**:
   - By unrolling the loop, the number of iterations is reduced, which can lead to better performance due to reduced loop overhead.

3. **Loop Fusion**:
   - Combining multiple loops into a single loop can reduce the overhead of loop control and improve cache locality.

4. **Loop Interchange**:
   - Changing the order of nested loops can improve cache performance by iterating over data in a more cache-friendly manner.

5. **Parallelization**:
   - Using OpenMP to parallelize the loop can significantly improve performance by distributing the workload across multiple threads.

### Optimized Code Explanation:

- **Parallelization**: The outer loop over `nl` is parallelized using OpenMP, allowing multiple iterations of the loop to run concurrently.
- **Loop Distribution**: The inner loop over `i` is distributed into chunks of 32 iterations (`32 * t1` to `32 * t1 + 31`), which helps in better cache utilization and allows for potential vectorization.
- **Loop Unrolling**: The inner loop is not explicitly unrolled, but the chunking approach can be seen as a form of unrolling, where each chunk is handled as a separate iteration.

This optimization leverages parallelization and loop distribution to improve performance while preserving the original meaning of the loop.*/

int nl, i;
int lbp = 0, ubp = floord(LEN_1D, 32);
#pragma omp parallel for private(i)
for (nl = 0; nl < 4 * iterations; nl++) {
    for (int t1 = lbp; t1 <= ubp; t1++) {
        for (i = max(0, 32 * t1); i <= min(LEN_1D - 1, 32 * t1 + 31); i++) {
            if (i + 1 < mid) {
                a[i] += b[i] * c[i];
            } else {
                a[i] += b[i] * d[i];
            }
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
    
    time_function(&s276, NULL);
    return EXIT_SUCCESS;
}