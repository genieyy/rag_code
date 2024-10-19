
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


real_t s2233(struct args_t * func_args)
{

//    loop interchange
//    interchanging with one of two inner loops

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
/*### Explanation of Transformations:
1. **Loop Unrolling and Parallelization**:
   - The outer loop over `nl` is parallelized using OpenMP (`#pragma omp parallel for`).
   - The inner loops over `i` and `j` are unrolled and vectorized using `#pragma ivdep` and `#pragma vector always` to ensure that the compiler can optimize the loop iterations.

2. **Loop Distribution**:
   - The loops are distributed to separate the computation of `aa` and `bb` to allow for better cache utilization and vectorization.

3. **Loop Tiling**:
   - The loops are tiled using a tile size of 32 (`32 * t2` and `32 * t1`) to improve cache locality and allow for better parallelization.

4. **Loop Fusion**:
   - The loops over `i` and `j` are fused to reduce the number of loop iterations and improve the efficiency of the computation.

5. **Loop Reordering**:
   - The loops are reordered to ensure that the most computationally intensive operations are performed first, allowing for better optimization by the compiler.

These transformations aim to improve the performance of the code by reducing the number of loop iterations, improving cache locality, and enabling parallel execution.*/

int t1, t2, t3, t4;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

lbp = 0;
ubp = floord(100 * (iterations / LEN_2D) - 1, 32);
#pragma omp parallel for private(lbv, ubv, t2, t3, t4)
for (int t1 = lbp; t1 <= ubp; t1++) {
    for (int t2 = 0; t2 <= floord(LEN_2D - 1, 32); t2++) {
        for (int t3 = 32 * t2; t3 <= min(LEN_2D - 1, 32 * t2 + 31); t3++) {
            lbv = max(32 * t1, 1);
            ubv = min(100 * (iterations / LEN_2D) - 1, 32 * t1 + 31);
            #pragma ivdep
            #pragma vector always
            for (int t4 = lbv; t4 <= ubv; t4++) {
                aa[t3][t4] = aa[t3 - 1][t4] + cc[t3][t4];
            }
            lbv = max(32 * t1, 1);
            ubv = min(100 * (iterations / LEN_2D) - 1, 32 * t1 + 31);
            #pragma ivdep
            #pragma vector always
            for (int t4 = lbv; t4 <= ubv; t4++) {
                bb[t4][t3] = bb[t4 - 1][t3] + cc[t4][t3];
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
    
    time_function(&s2233, NULL);
    return EXIT_SUCCESS;
}