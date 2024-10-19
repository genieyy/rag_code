
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


real_t s3251(struct args_t * func_args)
{

//    scalar and array expansion
//    scalar expansion

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
/*### Optimized Code



### Analysis of Loop Transformation Methods

1. **Loop Tiling/Blocking**:
   - The original loops are transformed into a nested structure where the outer loops iterate over blocks of the original loop indices. This is evident in the use of `floord` and `ceild` functions to determine the bounds of these blocks.
   - Example: `for (int t1 = lbp; t1 <= ubp; t1++)` where `lbp` and `ubp` are calculated to tile the loop over `LEN_1D`.

2. **Parallelization**:
   - The `#pragma omp parallel for` directive is used to parallelize the outer loop, allowing multiple threads to execute the loop iterations concurrently.
   - Example: `#pragma omp parallel for private(lbv, ubv)`

3. **Loop Fusion/Fission**:
   - The original loop body is split into multiple loops or fused into a single loop to improve cache locality and reduce loop overhead.
   - Example: The original loop body is kept within a single loop, but the outer loop is parallelized and tiled.

4. **Loop Unrolling**:
   - Although not explicitly shown in the provided examples, loop unrolling can be inferred as a potential optimization if the block size (`32` in this case) is chosen to align with the hardware's cache line size.

### Learning from the Examples

- **Parallelization**: Utilizing OpenMP for parallel execution can significantly improve performance by distributing the workload across multiple threads.
- **Loop Tiling**: Breaking the loop into smaller blocks can improve cache utilization, especially when dealing with large arrays.
- **Private Variables**: Ensuring that each thread has its own copy of loop variables (`lbv`, `ubv`) prevents race conditions and ensures thread safety.

### Application to the Given Code

- **Parallelization**: The outer loop over `iterations` is not parallelized because it involves a function call (`dummy`) that might have side effects. However, the inner loop over `i` is parallelized using OpenMP.
- **Loop Tiling**: The inner loop is tiled using `floord(LEN_1D - 2, 32)` to create blocks of size 32, which can help in better cache utilization.
- **Private Variables**: The loop variables `lbv` and `ubv` are declared as private to ensure thread safety.*/

int lb, ub, lbp, ubp;
register int lbv, ubv;
lbp = 0;
ubp = floord(LEN_1D - 2, 32);

for (int nl = 0; nl < iterations; nl++) {
#pragma omp parallel for private(lbv, ubv)
    for (int t1 = lbp; t1 <= ubp; t1++) {
        for (int t2 = max(0, 32 * t1); t2 <= min(LEN_1D - 2, 32 * t1 + 31); t2++) {
            a[t2 + 1] = b[t2] + c[t2];
            b[t2] = c[t2] * e[t2];
            d[t2] = a[t2] * e[t2];
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
    
    time_function(&s3251, NULL);
    return EXIT_SUCCESS;
}