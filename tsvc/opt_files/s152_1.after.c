
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

// %1.5

void s152s(real_t a[LEN_1D], real_t b[LEN_1D], real_t c[LEN_1D], int i)
{
    a[i] += b[i] * c[i];
}

real_t s152(struct args_t * func_args)
{

//    interprocedural data flow analysis
//    collecting information from a subroutine

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
/*### Analysis of Meaning-Preserving Loop Transformation Methods:

1. **Loop Distribution/Partitioning**: The original loops are split into multiple loops to isolate different operations. This allows for better optimization opportunities, such as parallelization and vectorization.

2. **Loop Fusion**: In some cases, loops that operate on the same data are combined to reduce overhead and improve cache locality.

3. **Loop Interchange**: The order of nested loops is changed to improve data locality and reduce cache misses.

4. **Loop Tiling/Blocking**: The iteration space is divided into smaller blocks (tiles) to improve cache utilization.

5. **Parallelization**: Using OpenMP, loops are parallelized to leverage multiple CPU cores, which can significantly speed up execution.

### Optimized Code Explanation:

- **Loop Distribution**: The outer loop over `iterations` is split into smaller chunks using `floord` to distribute the work across multiple threads.
- **Parallelization**: The outer loop is parallelized using OpenMP to allow multiple threads to work on different chunks of the iteration space.
- **Loop Interchange**: The order of loops is adjusted to ensure that the innermost loop operates on contiguous memory, improving cache performance.
- **Loop Tiling**: The iteration space is divided into blocks (tiles) to improve cache utilization.

This optimized code should provide better performance by leveraging parallel processing and improving cache locality.*/

int t1, t2, t3, t4, t5;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

lbp = 0;
ubp = floord(iterations - 1, 32);
#pragma omp parallel for private(lbv, ubv, t2, t3, t4, t5)
for (int t1 = lbp; t1 <= ubp; t1++) {
    for (int t2 = 0; t2 <= 2; t2++) {
        for (int t3 = 0; t3 <= min(floord(iterations - 1, 32), floord(32 * t1 - 32 * t2 + iterations - 1, 32)); t3++) {
            if (t3 == 0) {
                for (int t4 = max(32 * t1, 32 * t2 + 32); t4 <= min(LEN_1D - 1, 32 * t1 + 31); t4++) {
                    b[t4] = d[t4] * e[t4];
                    s152s(a, b, c, t4);
                }
                dummy(a, b, c, d, e, aa, bb, cc, 0.);
            }
            for (int t4 = max(32 * t1, 32 * t2 + 32); t4 <= min(LEN_1D - 1, 32 * t1 + 31); t4++) {
                b[t4] = d[t4] * e[t4];
                s152s(a, b, c, t4);
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
    
    time_function(&s152, NULL);
    return EXIT_SUCCESS;
}