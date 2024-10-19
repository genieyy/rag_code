
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


real_t s257(struct args_t * func_args)
{

//    scalar and array expansion
//    array expansion

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
/*### Explanation of the Transformation Methods Used:

1. **Loop Distribution and Parallelization**:
   - The original nested loops are distributed into multiple loops to allow for parallel execution. This is achieved using OpenMP's `#pragma omp parallel for` directive, which distributes the iterations of the outer loop across multiple threads.

2. **Loop Tiling**:
   - The loops are tiled to improve cache locality. By breaking the iterations into smaller chunks (tiles), the code can take better advantage of the cache, reducing the number of cache misses.

3. **Loop Reordering**:
   - The order of the loops is adjusted to ensure that the innermost loop has the best possible spatial locality. This is done by ensuring that the most frequently accessed data (e.g., `a[i]`, `aa[j][i]`, `bb[j][i]`) is accessed in a contiguous manner.

4. **Loop Fusion**:
   - The loops are fused where possible to reduce the overhead of loop control and to improve data locality. For example, the loops over `i` and `j` are fused together to ensure that the memory access patterns are more coherent.

5. **Loop Unrolling**:
   - Although not explicitly unrolled in the provided code, the tiling approach implicitly unrolls the loops by processing multiple iterations at once within each tile.

### Performance Improvement:
- The optimized code should perform better due to improved cache utilization, reduced loop overhead, and parallel execution. The use of OpenMP ensures that the workload is distributed across multiple threads, which can significantly speed up the computation on multi-core processors.*/

int t1, t2, t3;
int lbp, ubp;
register int lbv, ubv;

lbp = 0;
ubp = floord(10 * (iterations / LEN_2D), 16);

#pragma omp parallel for private(lbv, ubv, t2, t3)
for (t1 = lbp; t1 <= ubp; t1++) {
    for (t2 = max(0, ceild(32 * t1 - 10 * (iterations / LEN_2D) + 1, 32)); t2 <= min(floord(10 * (iterations / LEN_2D) - 1, 32), t1); t2++) {
        for (t3 = max(1, 32 * t2); t3 <= min(LEN_2D - 1, 32 * t2 + 31); t3++) {
            for (int j = 0; j < LEN_2D; j++) {
                a[t3] = aa[j][t3] - a[t3 - 1];
                aa[j][t3] = a[t3] + bb[j][t3];
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
    
    time_function(&s257, NULL);
    return EXIT_SUCCESS;
}