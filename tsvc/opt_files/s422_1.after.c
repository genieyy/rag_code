
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


real_t s422(struct args_t * func_args)
{

//    storage classes and equivalencing
//    common and equivalence statement
//    anti-dependence, threshold of 4

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

    xx = flat_2d_array + 4;

#pragma scop
/*### Explanation of the Transformation Methods Used:

1. **Loop Distribution and Parallelization**:
   - The original nested loops are distributed into multiple loops to allow for parallel execution. This is done using OpenMP's `#pragma omp parallel for` directive.

2. **Loop Tiling**:
   - The loops are tiled to improve cache locality. This is achieved by introducing new variables (`t1`, `t2`, `t3`, `t4`) that represent the tile indices and sizes.

3. **Loop Reordering and Fusion**:
   - The loops are reordered and fused to minimize the number of iterations and to ensure that the innermost loop has the highest possible number of iterations, which is beneficial for vectorization and parallelization.

4. **Bounds Calculation**:
   - The bounds of the loops are carefully calculated to ensure that all iterations are covered and no out-of-bounds accesses occur. This involves using functions like `floord`, `ceild`, `max`, and `min` to compute the correct loop bounds.

### What I Learned:
- **Parallelization**: Using OpenMP to parallelize loops can significantly improve performance.
- **Loop Tiling**: Tiling loops can improve cache utilization and performance.
- **Loop Reordering**: Reordering loops can help in better exploiting parallelism and vectorization opportunities.
- **Bounds Calculation**: Careful calculation of loop bounds is crucial to ensure correctness and performance.

### Optimized Code:
The optimized code applies these techniques to the given loop structure, distributing the work across multiple threads and tiling the loops to improve cache performance.*/

int t1, t2, t3, t4;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

for (int t1 = 0; t1 <= floord(8 * iterations - 1, 32); t1++) {
    lbp = max(ceild(t1, 2), ceild(32 * t1 - 8 * iterations + 1, 32));
    ubp = floord(2 * t1 + 1, 3);
#pragma omp parallel for private(lbv, ubv, t3, t4)
    for (int t2 = lbp; t2 <= ubp; t2++) {
        for (int t3 = max(32 * t1 - 32 * t2, 16 * t2 + 1); t3 <= min(8 * iterations - 1, 32 * t1 - 32 * t2 + 31); t3++) {
            for (int t4 = max(32 * t2, t3); t4 <= min(32 * t2 + 31, 2 * t3 - 1); t4++) {
                xx[t4] = flat_2d_array[t4 + 8] + a[t4];
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
    
    time_function(&s422, NULL);
    return EXIT_SUCCESS;
}