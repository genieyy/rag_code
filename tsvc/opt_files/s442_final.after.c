
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


real_t s442(struct args_t * func_args)
{

//    non-logical if's
//    computed goto

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
#include <immintrin.h>  // For AVX intrinsics

#define VECTOR_SIZE 4  // Assuming AVX2, which operates on 256-bit registers (4 doubles)

for (int nl = 0; nl < iterations/2; nl++) {
    for (int i = 0; i < LEN_1D; i += VECTOR_SIZE) {
        __m256d vec_a = _mm256_loadu_pd(&a[i]);
        __m256d vec_b = _mm256_loadu_pd(&b[i]);
        __m256d vec_c = _mm256_loadu_pd(&c[i]);
        __m256d vec_d = _mm256_loadu_pd(&d[i]);
        __m256d vec_e = _mm256_loadu_pd(&e[i]);
        __m256i vec_indx = _mm256_loadu_si256((__m256i*)&indx[i]);

        __m256d vec_b_sq = _mm256_mul_pd(vec_b, vec_b);
        __m256d vec_c_sq = _mm256_mul_pd(vec_c, vec_c);
        __m256d vec_d_sq = _mm256_mul_pd(vec_d, vec_d);
        __m256d vec_e_sq = _mm256_mul_pd(vec_e, vec_e);

        __m256d vec_mask1 = _mm256_castsi256_pd(_mm256_cmpeq_epi64(vec_indx, _mm256_set1_epi64x(1)));
        __m256d vec_mask2 = _mm256_castsi256_pd(_mm256_cmpeq_epi64(vec_indx, _mm256_set1_epi64x(2)));
        __m256d vec_mask3 = _mm256_castsi256_pd(_mm256_cmpeq_epi64(vec_indx, _mm256_set1_epi64x(3)));
        __m256d vec_mask4 = _mm256_castsi256_pd(_mm256_cmpeq_epi64(vec_indx, _mm256_set1_epi64x(4)));

        vec_a = _mm256_blendv_pd(vec_a, _mm256_add_pd(vec_a, vec_b_sq), vec_mask1);
        vec_a = _mm256_blendv_pd(vec_a, _mm256_add_pd(vec_a, vec_c_sq), vec_mask2);
        vec_a = _mm256_blendv_pd(vec_a, _mm256_add_pd(vec_a, vec_d_sq), vec_mask3);
        vec_a = _mm256_blendv_pd(vec_a, _mm256_add_pd(vec_a, vec_e_sq), vec_mask4);

        _mm256_storeu_pd(&a[i], vec_a);
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
    
    time_function(&s442, NULL);
    return EXIT_SUCCESS;
}