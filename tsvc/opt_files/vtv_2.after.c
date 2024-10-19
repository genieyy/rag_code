
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


real_t vtv(struct args_t * func_args)
{

//    control loops
//    vector times vector

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
/*The examples provided demonstrate several loop transformation techniques to optimize performance:

1. **Loop Tiling/Blocking**: This technique divides the iteration space into smaller blocks, which can help in better cache utilization and parallelization.
2. **Parallelization with OpenMP**: The use of `#pragma omp parallel for` directives to parallelize loops.
3. **Loop Fusion/Fission**: Combining or splitting loops to improve locality and reduce overhead.
4. **Loop Unrolling**: Reducing loop overhead by processing multiple elements per iteration.
5. **Loop Interchange**: Changing the order of nested loops to improve cache performance.

Given these techniques, we can apply similar optimizations to the provided loop:



### Explanation:
1. **Loop Tiling**: The iteration space of `LEN_1D` is divided into smaller blocks (tiles) of size `TILE_SIZE`. This helps in better cache utilization.
2. **Parallelization**: The outer loop over the tiles is parallelized using OpenMP, allowing multiple threads to process different tiles concurrently.
3. **Loop Fusion**: The inner loop over `i` is fused with the outer loop over `nl`, ensuring that the computation is done in a single loop structure.

This optimization should improve the performance of the original loop by leveraging parallelism and better cache utilization.*/

#include <omp.h>

#define TILE_SIZE 32

for (int nl = 0; nl < iterations * 10; nl++) {
    int len_1d_tiled = (LEN_1D + TILE_SIZE - 1) / TILE_SIZE;
    #pragma omp parallel for
    for (int tile = 0; tile < len_1d_tiled; tile++) {
        int start = tile * TILE_SIZE;
        int end = (tile + 1) * TILE_SIZE;
        if (end > LEN_1D) end = LEN_1D;
        for (int i = start; i < end; i++) {
            a[i] *= b[i];
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
    
    time_function(&vtv, NULL);
    return EXIT_SUCCESS;
}