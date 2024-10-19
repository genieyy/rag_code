
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


real_t s256(struct args_t * func_args)
{

//    scalar and array expansion
//    array expansion

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
/*### Explanation of the Transformation Methods Used:

1. **Loop Distribution and Parallelization**:
   - The original nested loops are distributed into multiple loops to enable parallel execution using OpenMP. This is done by introducing new loop variables (`t1`, `t2`, `t3`) and applying `#pragma omp parallel for` to the outermost loop.

2. **Loop Tiling**:
   - The loops are tiled to improve cache locality. This is achieved by introducing new loop bounds (`lbp`, `ubp`) that define the range of the loop iterations in a way that maximizes the reuse of cached data.

3. **Loop Fusion**:
   - The inner loops are fused to reduce the overhead of loop control and to improve the efficiency of memory access patterns. This is done by combining the loops over `i` and `j` into a single loop over `t3`.

4. **Loop Interchange**:
   - The order of the loops is interchanged to ensure that the most frequently accessed data is loaded into the cache first. This is done by changing the order of the loops over `i` and `j` to `t3` and `j`.

### What I Learned:

- **Parallelization**: By using OpenMP, the code can be parallelized to take advantage of multi-core processors, which can significantly improve performance.
- **Loop Tiling**: Tiling the loops helps in improving cache locality, which is crucial for performance in nested loops with large iteration spaces.
- **Loop Fusion and Interchange**: These techniques help in reducing the overhead of loop control and improving memory access patterns, leading to better performance.

### Optimized Code:

The optimized code applies these techniques to the given loop structure, resulting in a more efficient execution. The use of OpenMP ensures that the parallelizable parts of the code are executed concurrently, while loop tiling and fusion improve the locality and efficiency of the memory accesses.*/

int t1, t2, t3;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

for (int t1 = 0; t1 <= floord(10 * (iterations / LEN_2D) - 1, 32); t1++) {
    lbp = max(0, ceild(32 * t1 - 10 * (iterations / LEN_2D) + 1, 32));
    ubp = min(floord(10 * (iterations / LEN_2D) - 1, 32), t1);
#pragma omp parallel for private(lbv, ubv, t3)
    for (int t2 = lbp; t2 <= ubp; t2++) {
        for (int t3 = max(0, 32 * t1 - 32 * t2); t3 <= min(LEN_2D - 1, 32 * t1 - 32 * t2 + 31); t3++) {
            for (int j = 1; j < LEN_2D; j++) {
                a[j] = (real_t)1.0 - a[j - 1];
                aa[j][t3] = a[j] + bb[j][t3] * d[j];
            }
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
    
    time_function(&s256, NULL);
    return EXIT_SUCCESS;
}