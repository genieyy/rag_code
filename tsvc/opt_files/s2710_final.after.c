
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


real_t s2710(struct args_t * func_args)
{

//    control flow
//    scalar and vector ifs

    int x = *(int*)func_args->arg_info;

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
#include <immintrin.h> // For AVX intrinsics

for (int nl = 0; nl < iterations/2; nl++) {
    for (int i = 0; i < LEN_1D - 7; i += 8) {
        // Load data into AVX registers
        __m256 a_vec = _mm256_loadu_ps(&a[i]);
        __m256 b_vec = _mm256_loadu_ps(&b[i]);
        __m256 d_vec = _mm256_loadu_ps(&d[i]);
        __m256 e_vec = _mm256_loadu_ps(&e[i]);
        __m256 c_vec = _mm256_loadu_ps(&c[i]);

        // Compare a_vec and b_vec
        __m256 mask = _mm256_cmp_ps(a_vec, b_vec, _CMP_GT_OQ);

        // Compute a_vec += b_vec * d_vec where a_vec > b_vec
        __m256 a_update = _mm256_add_ps(a_vec, _mm256_mul_ps(b_vec, d_vec));
        a_vec = _mm256_blendv_ps(a_vec, a_update, mask);

        // Compute b_vec = a_vec + e_vec * e_vec where a_vec <= b_vec
        __m256 b_update = _mm256_add_ps(a_vec, _mm256_mul_ps(e_vec, e_vec));
        b_vec = _mm256_blendv_ps(b_update, b_vec, mask);

        // Compute c_vec based on conditions
        __m256 d_sq = _mm256_mul_ps(d_vec, d_vec);
        __m256 e_sq = _mm256_mul_ps(e_vec, e_vec);
        __m256 c_update1 = _mm256_add_ps(c_vec, d_sq);
        __m256 c_update2 = _mm256_add_ps(c_vec, e_sq);
        __m256 c_update3 = _mm256_add_ps(a_vec, d_sq);
        __m256 c_update4 = _mm256_add_ps(_mm256_set1_ps((real_t)1.), _mm256_mul_ps(d_vec, e_vec));

        // Update c_vec based on conditions
        __m256 c_update = _mm256_blendv_ps(c_update2, c_update1, mask);
        c_update = _mm256_blendv_ps(c_update4, c_update3, _mm256_cmp_ps(_mm256_set1_ps(x), _mm256_setzero_ps(), _CMP_GT_OQ));
        c_vec = _mm256_blendv_ps(c_vec, c_update, mask);

        // Store results back to memory
        _mm256_storeu_ps(&a[i], a_vec);
        _mm256_storeu_ps(&b[i], b_vec);
        _mm256_storeu_ps(&c[i], c_vec);
    }

    // Handle the remaining elements
    for (int i = LEN_1D - (LEN_1D % 8); i < LEN_1D; i++) {
        if (a[i] > b[i]) {
            a[i] += b[i] * d[i];
            if (LEN_1D > 10) {
                c[i] += d[i] * d[i];
            } else {
                c[i] = d[i] * e[i] + (real_t)1.;
            }
        } else {
            b[i] = a[i] + e[i] * e[i];
            if (x > (real_t)0.) {
                c[i] = a[i] + d[i] * d[i];
            } else {
                c[i] += e[i] * e[i];
            }
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
    
    time_function(&s2710, &s1);
    return EXIT_SUCCESS;
}