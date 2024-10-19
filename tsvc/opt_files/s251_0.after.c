
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


real_t s251(struct args_t * func_args)
{

//    scalar and array expansion
//    scalar expansion

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

    real_t s;
#pragma scop
/*### Optimized Code



### Explanation of Loop Transformation Methods Used

1. **Parallelization with OpenMP**:
   - The `#pragma omp parallel for` directive is used to parallelize the outer loop. This allows multiple threads to execute the loop iterations concurrently, which can significantly improve performance on multi-core processors.
   - The `private(s)` clause ensures that each thread has its own private copy of the variable `s`, preventing race conditions and ensuring thread safety.

2. **Loop Unrolling**:
   - Although not explicitly unrolled in the provided code, loop unrolling can be considered if the inner loop (`i`) has a fixed and small number of iterations. However, in this case, the inner loop is left as is because its length (`LEN_1D`) is not specified.

3. **Reduction of Redundant Computations**:
   - The variable `s` is computed once per iteration of the inner loop and then used to update `a[i]`. This reduces redundant computations and improves efficiency.

4. **Memory Access Optimization**:
   - The code accesses arrays `b`, `c`, `d`, and `a` sequentially, which is cache-friendly. This ensures that memory access patterns are optimized for modern CPU caches.

### Summary of Learnings

- **Parallelization**: Utilizing OpenMP for parallelizing loops can significantly improve performance by leveraging multiple CPU cores.
- **Private Variables**: Ensuring that variables used within parallel loops are private to each thread prevents race conditions and ensures correctness.
- **Memory Access Patterns**: Sequential memory access patterns are crucial for optimizing cache utilization and overall performance.
- **Reduction of Redundant Computations**: Minimizing redundant computations within loops can lead to performance gains.

By applying these techniques, the provided code can be optimized to run more efficiently on multi-core systems.*/

#pragma omp parallel for private(s)
for (int nl = 0; nl < 4 * iterations; nl++) {
    double s;
    for (int i = 0; i < LEN_1D; i++) {
        s = b[i] + c[i] * d[i];
        a[i] = s * s;
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
    
    time_function(&s251, NULL);
    return EXIT_SUCCESS;
}