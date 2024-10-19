
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


real_t s491(struct args_t * func_args)
{

//    vector semantics
//    indirect addressing on lhs, store in sequence
//    scatter is required

    int * __restrict__ ip = func_args->arg_info;

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
for (int nl = 0; nl < iterations; nl++) {
    // Temporary variables to store intermediate results
    double temp_a[4];
    double temp_b[4];
    double temp_c[4];
    double temp_d[4];
    int temp_ip[4];

    // Unroll the loop by 4 and use temporary variables
    for (int i = 0; i < LEN_1D - 3; i += 4) {
        temp_ip[0] = ip[i];
        temp_ip[1] = ip[i + 1];
        temp_ip[2] = ip[i + 2];
        temp_ip[3] = ip[i + 3];

        temp_b[0] = b[i];
        temp_b[1] = b[i + 1];
        temp_b[2] = b[i + 2];
        temp_b[3] = b[i + 3];

        temp_c[0] = c[i];
        temp_c[1] = c[i + 1];
        temp_c[2] = c[i + 2];
        temp_c[3] = c[i + 3];

        temp_d[0] = d[i];
        temp_d[1] = d[i + 1];
        temp_d[2] = d[i + 2];
        temp_d[3] = d[i + 3];

        temp_a[0] = temp_b[0] + temp_c[0] * temp_d[0];
        temp_a[1] = temp_b[1] + temp_c[1] * temp_d[1];
        temp_a[2] = temp_b[2] + temp_c[2] * temp_d[2];
        temp_a[3] = temp_b[3] + temp_c[3] * temp_d[3];

        a[temp_ip[0]] = temp_a[0];
        a[temp_ip[1]] = temp_a[1];
        a[temp_ip[2]] = temp_a[2];
        a[temp_ip[3]] = temp_a[3];
    }

    // Handle the remaining elements
    for (int i = LEN_1D - (LEN_1D % 4); i < LEN_1D; i++) {
        a[ip[i]] = b[i] + c[i] * d[i];
    }

    // Call the dummy function once per iteration
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
    
    time_function(&s491, ip);
    return EXIT_SUCCESS;
}