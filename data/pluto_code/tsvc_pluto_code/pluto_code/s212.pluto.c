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



// array definitions
__attribute__((aligned(ARRAY_ALIGNMENT))) real_t flat_2d_array[LEN_2D*LEN_2D];

__attribute__((aligned(ARRAY_ALIGNMENT))) real_t x[LEN_1D];

__attribute__((aligned(ARRAY_ALIGNMENT))) real_t a[LEN_1D],b[LEN_1D],c[LEN_1D],d[LEN_1D],e[LEN_1D],
                                   aa[LEN_2D][LEN_2D],bb[LEN_2D][LEN_2D],cc[LEN_2D][LEN_2D],tt[LEN_2D][LEN_2D];

__attribute__((aligned(ARRAY_ALIGNMENT))) int indx[LEN_1D];

real_t* __restrict__ xx;
real_t* yy;


real_t s212(struct args_t * func_args)
{

//    statement reordering
//    dependency needing temporary

    initialise_arrays(__func__);
    func_args->c1 = omp_get_wtime();
  int t1, t2, t3, t4, t5;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
for (t1=0;t1<=floord(2*iterations+LEN_1D-4,32);t1++) {
  lbp=max(ceild(t1,2),ceild(32*t1-iterations+1,32));
  ubp=min(floord(32*t1+LEN_1D+29,64),t1);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5)
  for (t2=lbp;t2<=ubp;t2++) {
    if (t1 == 2*t2) {
      for (t3=16*t1;t3<=min(iterations-1,16*t1+31);t3++) {
        if (t1%2 == 0) {
          b[0] += a[0 + 1] * d[0];;
        }
        if (t1%2 == 0) {
          a[0] *= c[0];;
        }
        if (t1%2 == 0) {
          dummy(a, b, c, d, e, aa, bb, cc, 0.);;
        }
        for (t4=t3+1;t4<=16*t1+31;t4++) {
          if (t1%2 == 0) {
            b[(-t3+t4)] += a[(-t3+t4) + 1] * d[(-t3+t4)];;
          }
          if (t1%2 == 0) {
            a[(-t3+t4)] *= c[(-t3+t4)];;
          }
        }
      }
    }
    if (t1 <= 2*t2-1) {
      for (t3=max(32*t1-32*t2,32*t2-LEN_1D+2);t3<=min(iterations-1,32*t1-32*t2+31);t3++) {
        for (t4=32*t2;t4<=min(32*t2+31,t3+LEN_1D-2);t4++) {
          b[(-t3+t4)] += a[(-t3+t4) + 1] * d[(-t3+t4)];;
          a[(-t3+t4)] *= c[(-t3+t4)];;
        }
      }
    }
  }
}
/* End of CLooG code */

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
}

int main(int argc, char ** argv){
    int n1 = 1;
    int n3 = 1;
    int* ip;
    real_t s1,s2;
    init(&ip, &s1, &s2);
    printf("Loop \tTime(sec) \tChecksum\n");
    
    time_function(&s212, NULL);
    return EXIT_SUCCESS;
}
