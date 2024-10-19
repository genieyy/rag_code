
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


real_t s4114(struct args_t * func_args)
{

//    indirect addressing
//    mix indirect addressing with variable lower and upper bounds
//    gather is required

    struct{int * __restrict__ a;int b;} * x = func_args->arg_info;
    int * __restrict__ ip = x->a;
    int n1 = x->b;

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

    int k;
#pragma scop
/*### Explanation:
- **Loop Unrolling**: The loop is unrolled to process 4 elements at a time, reducing the overhead of loop control.
- **Reduced Indirection**: The indices `k1`, `k2`, `k3`, and `k4` are precomputed for each iteration of the unrolled loop, reducing the number of indirection operations. This can help in improving cache locality and reducing the number of memory accesses.
- **Remaining Elements**: The loop handles the remaining elements that do not fit into the unrolled block.*/

/* Optimized Code 2: Loop Unrolling with Reduced Indirection
This version combines loop unrolling with reducing the number of indirection operations by precomputing the indices.

*/
for (int nl = 0; nl < iterations; nl++) {
    int i = n1-1;
    for (; i + 3 < LEN_1D; i += 4) {
        int k1 = ip[i];
        int k2 = ip[i+1];
        int k3 = ip[i+2];
        int k4 = ip[i+3];

        a[i] = b[i] + c[LEN_1D-k1+1-2] * d[i];
        a[i+1] = b[i+1] + c[LEN_1D-k2+1-2] * d[i+1];
        a[i+2] = b[i+2] + c[LEN_1D-k3+1-2] * d[i+2];
        a[i+3] = b[i+3] + c[LEN_1D-k4+1-2] * d[i+3];

        k1 += 5;
        k2 += 5;
        k3 += 5;
        k4 += 5;
    }
    for (; i < LEN_1D; i++) {
        k = ip[i];
        a[i] = b[i] + c[LEN_1D-k+1-2] * d[i];
        k += 5;
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
    
    time_function(&s4114, &(struct{int*a;int b;}){ip, n1});
    return EXIT_SUCCESS;
}