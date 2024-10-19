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


real_t s113(struct args_t * func_args)
{

//    linear dependence testing
//    a(i)=a(1) but no actual dependence cycle

    initialise_arrays(__func__);
    func_args->c1 = omp_get_wtime();
  int t1, t2, t3, t4;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
lbp=0;
ubp=floord(4*iterations+LEN_1D-32004,32);
#pragma omp parallel for private(lbv,ubv,t2,t3,t4)
for (t1=lbp;t1<=ubp;t1++) {
  for (t2=0;t2<=min(min(floord(iterations-1,8),floord(-799992*t1+99999*iterations-99999,735992)),floord(-3200000*t1+368001*iterations+100000*LEN_1D-468001,2944008));t2++) {
    if (t2 >= 1) {
      for (t3=max(1,32*t1);t3<=min(LEN_1D-1,32*t1+31);t3++) {
        for (t4=32*t2;t4<=min(4*iterations-1,32*t2+31);t4++) {
          a[t3] = a[0] + b[t3];;
        }
      }
    }
    if ((t1 == 0) && (t2 == 0)) {
      dummy(a, b, c, d, e, aa, bb, cc, 0.);;
    }
    if (t2 == 0) {
      for (t3=max(1,32*t1);t3<=min(LEN_1D-1,32*t1+31);t3++) {
        a[t3] = a[0] + b[t3];;
        dummy(a, b, c, d, e, aa, bb, cc, 0.);;
        for (t4=1;t4<=31;t4++) {
          a[t3] = a[0] + b[t3];;
        }
      }
    }
    if (t2 == 0) {
      for (t3=max(LEN_1D,32*t1);t3<=min(4*iterations-1,32*t1+31);t3++) {
        dummy(a, b, c, d, e, aa, bb, cc, 0.);;
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
    
    time_function(&s113, NULL);
    return EXIT_SUCCESS;
}
