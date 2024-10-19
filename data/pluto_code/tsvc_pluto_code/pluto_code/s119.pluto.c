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


real_t s119(struct args_t * func_args)
{

//    linear dependence testing
//    no dependence - vectorizable

    initialise_arrays(__func__);
    func_args->c1 = omp_get_wtime();
  int t1, t2, t3, t4, t5, t6;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
for (t1=0;t1<=floord(25*iterations-1,4);t1++) {
  for (t2=0;t2<=floord(32*t1+LEN_2D+30,32);t2++) {
    for (t3=max(max(max(max(max(max(0,ceild(10000227*t1+10000096*t2-62500600*iterations-312503*LEN_2D+625006,10000227)),ceild(143*t1+10000096*t2-89375858,10000239)),ceild(20000455*t1+20000192*t2-125001200*iterations-159376530,20000455)),ceild(285*t1+20000192*t2-625006*LEN_2D-18125174,20000477)),ceild(380*t1+26666656*t2-2375*iterations-833333*LEN_2D+213333343,26667036)),ceild(572*t1+39999984*t2-3575*iterations+143,40000556));t3<=min(min(min(min(min(min(floord(32*t1+LEN_2D+30,32),floord(-143*t1+10000239*t2+89375858,10000096)),floord(-10000227*t1+10000227*t2+62500600*iterations+312503*LEN_2D-625006,10000096)),floord(-20000455*t1+20000455*t2+125001200*iterations+159376530,20000192)),floord(-285*t1+20000477*t2+625006*LEN_2D+18125174,20000192)),floord(-380*t1+26667036*t2+2375*iterations+833333*LEN_2D-213333343,26666656)),floord(-572*t1+40000556*t2+3575*iterations-143,39999984));t3++) {
      if ((t1 == 0) && (t2 == 0) && (t3 == 0)) {
        for (t4=0;t4<=30;t4++) {
          dummy(a, b, c, d, e, aa, bb, cc, 0.);;
          for (t5=t4+1;t5<=31;t5++) {
            lbp=t4+1;
            ubp=31;
#pragma omp parallel for private(lbv,ubv)
            for (t6=lbp;t6<=ubp;t6++) {
              aa[(-t4+t5)][(-t4+t6)] = aa[(-t4+t5)-1][(-t4+t6)-1] + bb[(-t4+t5)][(-t4+t6)];;
            }
          }
        }
      }
      if ((t1 == 0) && (t2 >= 1) && (t3 == 0)) {
        for (t4=max(0,32*t2-LEN_2D+1);t4<=30;t4++) {
          for (t5=32*t2;t5<=min(32*t2+31,t4+LEN_2D-1);t5++) {
            lbp=t4+1;
            ubp=31;
#pragma omp parallel for private(lbv,ubv)
            for (t6=lbp;t6<=ubp;t6++) {
              aa[(-t4+t5)][(-t4+t6)] = aa[(-t4+t5)-1][(-t4+t6)-1] + bb[(-t4+t5)][(-t4+t6)];;
            }
          }
        }
      }
      if (t3 >= 1) {
        for (t4=max(max(32*t1,32*t2-LEN_2D+1),32*t3-LEN_2D+1);t4<=min(min(min(200*iterations-1,32*t1+31),32*t2+30),32*t3+30);t4++) {
          for (t5=max(32*t2,t4+1);t5<=min(32*t2+31,t4+LEN_2D-1);t5++) {
            lbp=max(32*t3,t4+1);
            ubp=min(32*t3+31,t4+LEN_2D-1);
#pragma omp parallel for private(lbv,ubv)
            for (t6=lbp;t6<=ubp;t6++) {
              aa[(-t4+t5)][(-t4+t6)] = aa[(-t4+t5)-1][(-t4+t6)-1] + bb[(-t4+t5)][(-t4+t6)];;
            }
          }
        }
      }
      if ((t2 == 0) && (t3 == 0)) {
        for (t4=max(31,32*t1);t4<=min(200*iterations-1,32*t1+31);t4++) {
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
    
    time_function(&s119, NULL);
    return EXIT_SUCCESS;
}
