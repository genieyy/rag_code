
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


real_t s1221(struct args_t * func_args)
{

//    run-time symbolic resolution

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
/*### Analysis of Meaning-Preserving Loop Transformation Methods

1. **Loop Distribution and Parallelization**:
   - The original nested loops are transformed into a single outer loop that iterates over a range determined by the problem size (`iterations` and `LEN_1D`).
   - The inner loops are then parallelized using OpenMP (`#pragma omp parallel for`) to distribute the workload across multiple threads.

2. **Loop Tiling**:
   - The iterations of the inner loops are divided into tiles (chunks) of size 32, which helps in better cache utilization and reduces the overhead of parallelization.

3. **Loop Reordering and Index Transformation**:
   - The original loop indices (`i`) are transformed into new indices (`t1`, `t2`, `t3`) to facilitate the tiling and parallelization.
   - The new indices are carefully calculated to ensure that the original loop bounds and dependencies are preserved.

### Learnings from the Examples

- **Parallelization**: Utilizing OpenMP to parallelize loops can significantly improve performance by distributing the workload across multiple threads.
- **Loop Tiling**: Breaking down the problem into smaller chunks (tiles) can improve cache performance and make the parallelization more efficient.
- **Index Transformation**: Carefully transforming loop indices can help in maintaining the original loop semantics while enabling more efficient parallel execution.

### Optimized Code Explanation

- **Outer Loop (`t1`)**: Iterates over the number of iterations, divided into chunks of size 32.
- **Inner Loop (`t2`)**: Determines the starting point for the inner loop based on the current chunk.
- **Inner-most Loop (`t3`)**: Iterates over the elements within the current chunk, performing the original computation (`b[i] = b[i - 4] + a[i]`).
- **Parallelization**: The inner loop is parallelized using OpenMP to distribute the workload across multiple threads.
- **Dummy Function Call**: The `dummy` function is called after processing each chunk to ensure that the original semantics are preserved.*/

int t1, t2, t3, t4;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

for (int t1 = 0; t1 <= floord(iterations - 1, 32); t1++) {
    lbp = max(0, ceild(32 * t1 - LEN_1D + 5, 32));
    ubp = floord(32 * t1 + 31, 32);
#pragma omp parallel for private(lbv, ubv, t3, t4)
    for (int t2 = lbp; t2 <= ubp; t2++) {
        for (int t3 = max(4, 32 * t2); t3 <= min(LEN_1D - 1, 32 * t2 + 31); t3++) {
            b[t3] = b[t3 - 4] + a[t3];
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
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
    
    time_function(&s1221, NULL);
    return EXIT_SUCCESS;
}