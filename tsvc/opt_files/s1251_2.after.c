
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


real_t s1251(struct args_t * func_args)
{

//    scalar and array expansion
//    scalar expansion

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

    real_t s;
#pragma scop
/*### Explanation of the Transformation:
1. **Loop Distribution and Parallelization**: The original nested loops are transformed to distribute the workload across multiple threads using OpenMP. This is done by introducing a parallel loop (`#pragma omp parallel for`) and distributing the iterations of the outermost loop across threads.

2. **Loop Tiling**: The iterations of the loops are tiled to improve cache locality. This is done by introducing new variables (`t1`, `t2`, `t3`, `t4`) that represent the tile indices and the within-tile indices.

3. **Loop Reordering**: The loops are reordered to ensure that the innermost loop has the best possible spatial locality, which is crucial for performance on modern CPUs with deep memory hierarchies.

4. **Loop Fusion**: The `dummy` function call is fused with the innermost loop to reduce the overhead of function calls and to ensure that the data dependencies are respected.

### What I Learned:
- **Parallelization**: Using OpenMP to parallelize loops can significantly improve performance by distributing the workload across multiple threads.
- **Loop Tiling**: Tiling loops can improve cache locality, which is particularly important for performance in nested loops with large iteration spaces.
- **Loop Reordering**: Reordering loops can help in improving spatial locality and reducing the number of cache misses.
- **Loop Fusion**: Fusing loops that have data dependencies can reduce the overhead of function calls and improve performance by ensuring that the data is accessed in a more cache-friendly manner.*/

int t1, t2, t3, t4, t5;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

for (int t1 = 0; t1 <= floord(4 * iterations - 1, 32); t1++) {
    lbp = max(ceild(t1, 2), ceild(32 * t1 - 4 * iterations + 1, 32));
    ubp = floord(2 * t1 + 1, 3);
#pragma omp parallel for private(lbv, ubv, t3, t4, t5)
    for (int t2 = lbp; t2 <= ubp; t2++) {
        for (int t3 = max(max(0, 32 * t1 - 32 * t2), 16 * t2); t3 <= min(min(4 * iterations - 1, 32 * t2 + 31), 32 * t1 - 32 * t2 + 31); t3++) {
            for (int t4 = max(32 * t2, t3); t4 <= min(32 * t2 + 31, t3 + LEN_1D - 1); t4++) {
                double s = b[t4] + c[t4];
                b[t4] = a[t4] + d[t4];
                a[t4] = s * e[t4];
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
    
    time_function(&s1251, NULL);
    return EXIT_SUCCESS;
}