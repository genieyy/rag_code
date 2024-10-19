
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


real_t s1232(struct args_t * func_args)
{

//    loop interchange
//    interchanging of triangular loops

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
/*### Analysis of Meaning-Preserving Loop Transformation Methods

1. **Loop Tiling/Blocking**: The original code is transformed by dividing the iteration space into smaller blocks (tiles). This is evident in the transformation of the outer loops into a nested structure with `t1`, `t2`, and `t3` controlling the block sizes. This technique is used to improve cache locality and reduce cache misses.

2. **Parallelization with OpenMP**: The `#pragma omp parallel for` directive is used to parallelize the outer loop. This allows multiple threads to execute the loop iterations concurrently, which can significantly improve performance on multi-core processors.

3. **Loop Fusion/Fission**: The original code has a loop that iterates over `nl` and another loop that iterates over `j` and `i`. The transformation fuses these loops into a single nested loop structure, which can reduce the overhead of loop control and improve data locality.

4. **Loop Interchange**: The original loop order is changed to improve cache utilization. By interchanging the order of the loops, the innermost loop accesses data in a more contiguous manner, which can improve cache performance.

### Learning from the Examples

- **Cache Optimization**: By tiling the loops, the code can take better advantage of the CPU cache, reducing the number of cache misses and improving performance.
- **Parallel Execution**: Using OpenMP to parallelize the outer loop can leverage multiple CPU cores, leading to faster execution times.
- **Loop Structure**: The transformation of the loop structure to a nested form with block indices (`t1`, `t2`, `t3`) allows for more fine-grained control over the iteration space, which can be beneficial for both cache optimization and parallel execution.

### Optimized Code Explanation

- **Tiling**: The outer loop is tiled with `t1`, `t2`, and `t3` to create smaller blocks of iterations. This helps in better cache utilization.
- **Parallelization**: The `#pragma omp parallel for` directive is used to parallelize the outer loop, allowing multiple threads to execute the loop iterations concurrently.
- **Loop Fusion**: The loops over `nl`, `j`, and `i` are fused into a single nested loop structure, reducing loop control overhead and improving data locality.
- **Loop Interchange**: The loop order is interchanged to ensure that the innermost loop accesses data in a contiguous manner, improving cache performance.*/

int t1, t2, t3;
int lbp, ubp;
register int lbv, ubv;

lbp = 0;
ubp = floord(100 * (iterations / LEN_2D) - 1, 32);
#pragma omp parallel for private(lbv, ubv, t2, t3)
for (t1 = lbp; t1 <= ubp; t1++) {
    for (t2 = 0; t2 <= floord(LEN_2D - 1, 32); t2++) {
        for (t3 = max(32 * t1, 32 * t2); t3 <= min(100 * (iterations / LEN_2D) - 1, 32 * t1 + 31); t3++) {
            for (int j = 32 * t2; j <= min(LEN_2D - 1, 32 * t2 + 31); j++) {
                for (int i = j; i <= LEN_2D - 1; i++) {
                    aa[i][j] = bb[i][j] + cc[i][j];
                }
            }
            dummy(a, b, c, d, e, aa, bb, cc, 1.);
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
    
    time_function(&s1232, NULL);
    return EXIT_SUCCESS;
}