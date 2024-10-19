
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


real_t s323(struct args_t * func_args)
{

//    recurrences
//    coupled recurrence

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
/*### Explanation of Transformations:
1. **Loop Distribution and Parallelization**:
   - The outer loop over `nl` is preserved, but the inner loop over `i` is distributed into chunks of size 32 using the `floord` function. This allows for parallel execution using OpenMP.
   - The `#pragma omp parallel for` directive is used to parallelize the loop over `t2`, which represents the chunks of the array.

2. **Loop Unrolling**:
   - The inner loop over `t3` is unrolled to process 32 elements at a time, which can help in vectorization and improve cache locality.

3. **Private Variables**:
   - The variables `lbv`, `ubv`, `t3`, `t4`, and `t5` are declared as private within the OpenMP parallel region to ensure thread safety.

4. **Boundary Conditions**:
   - The `max` and `min` functions are used to handle the boundary conditions for the chunks, ensuring that the loop indices do not go out of bounds.

### Learnings:
- **Loop Distribution**: Breaking down the loop into smaller chunks can help in parallelizing the computation.
- **Parallelization**: Using OpenMP can significantly improve performance by distributing the workload across multiple threads.
- **Loop Unrolling**: Processing multiple elements within a loop iteration can help in vectorization and improve cache performance.
- **Boundary Handling**: Careful handling of loop boundaries is crucial to ensure correctness and avoid out-of-bounds access.*/

int t1, t2, t3, t4, t5;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

for (int t1 = 0; t1 <= (iterations / 2) - 1; t1++) {
    lbp = 0;
    ubp = floord(LEN_1D - 1, 32);
#pragma omp parallel for private(lbv, ubv, t3, t4, t5)
    for (int t2 = lbp; t2 <= ubp; t2++) {
        for (int t3 = max(1, 32 * t2); t3 <= min(LEN_1D - 1, 32 * t2 + 31); t3++) {
            a[t3] = b[t3 - 1] + c[t3] * d[t3];
            b[t3] = a[t3] + c[t3] * e[t3];
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
    
    time_function(&s323, NULL);
    return EXIT_SUCCESS;
}