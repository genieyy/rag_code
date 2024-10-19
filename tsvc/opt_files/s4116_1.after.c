
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


real_t s4116(struct args_t * func_args)
{

//    indirect addressing
//    more complicated sparse sdot
//    gather is required

    struct{int * __restrict__ a;int b;int c;} * x = func_args->arg_info;
    int * __restrict__ ip = x->a;
    int j = x->b;
    int inc = x->c;

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

    real_t sum;
    int off;
#pragma scop
```c
/* Optimized Code 1: Unroll Inner Loop
Unroll the inner loop to reduce loop overhead and improve instruction-level parallelism.
*/
for (int nl = 0; nl < 100 * iterations; nl++) {
    sum = 0.;
    for (int i = 0; i < LEN_2D - 1; i += 2) {
        off = inc + i;
        sum += a[off] * aa[j-1][ip[i]];
        off = inc + i + 1;
        sum += a[off] * aa[j-1][ip[i + 1]];
    }
    dummy(a, b, c, d, e, aa, bb, cc, 0.);
}
```

```c
/* Optimized Code 2: Reduce Redundant Calculations
Reduce redundant calculations by precomputing values that are used multiple times within the loop.
*/
for (int nl = 0; nl < 100 * iterations; nl++) {
    sum = 0.;
    double *aa_j_1 = aa[j-1];
    for (int i = 0; i < LEN_2D - 1; i++) {
        off = inc + i;
        sum += a[off] * aa_j_1[ip[i]];
    }
    dummy(a, b, c, d, e, aa, bb, cc, 0.);
}
```

```c
/* Optimized Code 3: Loop Fusion
Fuse the outer and inner loops to reduce loop overhead and improve cache locality.
*/
for (int nl = 0; nl < 100 * iterations * (LEN_2D - 1); nl++) {
    int i = nl % (LEN_2D - 1);
    int nl_outer = nl / (LEN_2D - 1);
    sum = 0.;
    off = inc + i;
    sum += a[off] * aa[j-1][ip[i]];
    if (i == LEN_2D - 2) {
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }
}
```#pragma endscop

    func_args->c2 = omp_get_wtime();
    return sum;
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
    
    time_function(&s4116, &(struct{int * a; int b; int c;}){ip, LEN_2D/2, n1});
    return EXIT_SUCCESS;
}