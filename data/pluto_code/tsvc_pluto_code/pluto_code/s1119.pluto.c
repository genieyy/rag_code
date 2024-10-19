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


real_t s1119(struct args_t * func_args)
{

//    linear dependence testing
//    no dependence - vectorizable

    initialise_arrays(__func__);
    func_args->c1 = omp_get_wtime();
  int t1, t2, t3, t4, t5, t6;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
lbp=0;
ubp=floord(200*iterations+LEN_2D-264,32);
#pragma omp parallel for private(lbv,ubv,t2,t3,t4,t5,t6)
for (t1=lbp;t1<=ubp;t1++) {
  for (t2=0;t2<=min(min(min(min(floord(25*iterations-1,4),floord(-1250012*t1+7812475*iterations+9687597,1249996)),floord(-9999996*t1+62499975*iterations-2499999,9999868)),floord(-26666656*t1+166666600*iterations+833333*LEN_2D-219999912,26666316)),floord(-80000768*t1+499998425*iterations+2500024*LEN_2D-22499961,79999748));t2++) {
    for (t3=t2;t3<=min(min(min(min(min(min(min(floord(-640008160*t1+4000051000*iterations+20000255*LEN_2D-5280067320,639991584),floord(-640008128*t1+4000050800*iterations+39999991*LEN_2D-10399999728,639991584)),floord(-9184*t1+639991584*t2+57400*iterations+287*LEN_2D-75768,639991584)),floord(-9152*t1+639991584*t2+57200*iterations+20000023*LEN_2D-5120008176,639991584)),floord(-640014560*t1+3999987400*iterations+20000455*LEN_2D+5079983480,639997984)),floord(-640014528*t1+3999987400*iterations+40000391*LEN_2D-60000328,639997984)),floord(-9184*t1+639997984*t2+287*LEN_2D+5739981632,639997984)),floord(-9152*t1+639997984*t2+20000223*LEN_2D+599997824,639997984));t3++) {
      if ((t2 == 0) && (t3 == 0)) {
        for (t4=32*t1;t4<=min(LEN_2D-1,32*t1+31);t4++) {
          dummy(a, b, c, d, e, aa, bb, cc, 0.);;
          for (t5=0;t5<=30;t5++) {
            for (t6=t5+1;t6<=31;t6++) {
              aa[(-t5+t6)][t4] = aa[(-t5+t6)-1][t4] + bb[(-t5+t6)][t4];;
            }
          }
        }
      }
      if (t3 >= 1) {
        for (t4=32*t1;t4<=min(LEN_2D-1,32*t1+31);t4++) {
          for (t5=max(32*t2,32*t3-LEN_2D+1);t5<=min(min(200*iterations-1,32*t2+31),32*t3+30);t5++) {
            for (t6=max(32*t3,t5+1);t6<=min(32*t3+31,t5+LEN_2D-1);t6++) {
              aa[(-t5+t6)][t4] = aa[(-t5+t6)-1][t4] + bb[(-t5+t6)][t4];;
            }
          }
        }
      }
      if ((t2 == 0) && (t3 == 0)) {
        for (t4=max(LEN_2D,32*t1);t4<=min(200*iterations-1,32*t1+31);t4++) {
          dummy(a, b, c, d, e, aa, bb, cc, 0.);;
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
    
    time_function(&s1119, NULL);
    return EXIT_SUCCESS;
}
