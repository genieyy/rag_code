
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
/*### Explanation of Optimizations:
1. **Loop Unrolling and Vectorization**: The inner loop is unrolled by a factor of 32, which helps in vectorization and reduces the loop overhead. This is similar to the loop transformations seen in the provided examples.
2. **Padding for Alignment**: Arrays `a`, `b`, `c`, and `d` are padded to ensure they are multiples of 32 elements long, which helps in better memory alignment and vectorization.
3. **Parallelization**: The outer loop is parallelized using OpenMP to leverage multiple cores, similar to the parallelization strategy in the provided examples.
4. **Conditional Check Reduction**: The conditional check `if (idx < LEN_1D)` is used to ensure that we do not access out-of-bounds memory in the padded region. This reduces the number of conditional checks inside the inner loop.*/

#include <omp.h>

#define LEN_1D_PADDED (LEN_1D + (LEN_1D % 32 == 0 ? 0 : 32 - (LEN_1D % 32)))

double a_padded[LEN_1D_PADDED];
double b_padded[LEN_1D_PADDED];
double c_padded[LEN_1D_PADDED];
double d_padded[LEN_1D_PADDED];

for (int nl = 0; nl < 4 * iterations; nl++) {
    #pragma omp parallel for
    for (int i = 0; i < LEN_1D_PADDED; i += 32) {
        for (int j = 0; j < 32; j++) {
            int idx = i + j;
            if (idx < LEN_1D) {
                if (idx + 1 < mid) {
                    a_padded[idx] += b_padded[idx] * c_padded[idx];
                } else {
                    a_padded[idx] += b_padded[idx] * d_padded[idx];
                }
            }
        }
    }
    dummy(a_padded, b_padded, c_padded, d_padded, e, aa, bb, cc, 0.);
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