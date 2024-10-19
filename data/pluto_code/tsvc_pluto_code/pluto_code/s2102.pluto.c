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


real_t s2102(struct args_t * func_args)
{

//    diagonals
//    identity matrix, best results vectorize both inner and outer loops

    initialise_arrays(__func__);
    func_args->c1 = omp_get_wtime();
  int t1, t2, t3, t4, t5, t6, t7;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
lbp=0;
ubp=floord(100*iterations+LEN_2D-260,32);
#pragma omp parallel for private(lbv,ubv,t2,t3,t4,t5,t6,t7)
for (t1=lbp;t1<=ubp;t1++) {
  for (t2=0;t2<=min(min(min(min(floord(25*iterations-1,8),floord(-2500024*t1+7812375*iterations+19687697,2499960)),floord(-19999992*t1+62499975*iterations-2499999,19999480)),floord(-26666656*t1+83333300*iterations+833333*LEN_2D-216666580,26665976)),floord(-80000768*t1+249996025*iterations+2500024*LEN_2D-12499865,79998728));t2++) {
    for (t3=0;t3<=min(min(min(floord(-256*t1+800*iterations+8*LEN_2D-2080,9999741),floord(-256*t1+8*LEN_2D+79998720,9999841)),floord(-255*t1+312503*LEN_2D-312503,9999841)),floord(-680*t1+2125*iterations+833333*LEN_2D-213333333,26665976));t3++) {
      if ((t1 == 0) && (t2 == 0) && (t3 == 0)) {
        for (t4=0;t4<=31;t4++) {
          dummy(a, b, c, d, e, aa, bb, cc, 0.);;
          for (t5=0;t5<=31;t5++) {
            lbv=0;
            ubv=t4-1;
#pragma ivdep
#pragma vector always
            for (t6=lbv;t6<=ubv;t6++) {
              aa[t6][t4] = 0.;;
            }
            aa[t4][t4] = 0.;;
            aa[t4][t4] = 1.;;
            lbv=t4+1;
            ubv=31;
#pragma ivdep
#pragma vector always
            for (t6=lbv;t6<=ubv;t6++) {
              aa[t6][t4] = 0.;;
            }
          }
        }
      }
      if ((t1 == 0) && (t2 >= 1) && (t3 == 0)) {
        for (t4=0;t4<=31;t4++) {
          for (t5=32*t2;t5<=min(100*iterations-1,32*t2+31);t5++) {
            lbv=0;
            ubv=t4-1;
#pragma ivdep
#pragma vector always
            for (t6=lbv;t6<=ubv;t6++) {
              aa[t6][t4] = 0.;;
            }
            aa[t4][t4] = 0.;;
            aa[t4][t4] = 1.;;
            lbv=t4+1;
            ubv=31;
#pragma ivdep
#pragma vector always
            for (t6=lbv;t6<=ubv;t6++) {
              aa[t6][t4] = 0.;;
            }
          }
        }
      }
      if ((t1 >= 1) && (t1 == t3)) {
        for (t4=32*t1;t4<=min(LEN_2D-1,32*t1+31);t4++) {
          for (t5=32*t2;t5<=min(100*iterations-1,32*t2+31);t5++) {
            lbv=32*t1;
            ubv=t4-1;
#pragma ivdep
#pragma vector always
            for (t6=lbv;t6<=ubv;t6++) {
              aa[t6][t4] = 0.;;
            }
            aa[t4][t4] = 0.;;
            aa[t4][t4] = 1.;;
            lbv=t4+1;
            ubv=min(LEN_2D-1,32*t1+31);
#pragma ivdep
#pragma vector always
            for (t6=lbv;t6<=ubv;t6++) {
              aa[t6][t4] = 0.;;
            }
          }
        }
      }
      if (t1 <= t3-1) {
        for (t4=32*t1;t4<=32*t1+31;t4++) {
          for (t5=32*t2;t5<=min(100*iterations-1,32*t2+31);t5++) {
            lbv=32*t3;
            ubv=min(LEN_2D-1,32*t3+31);
#pragma ivdep
#pragma vector always
            for (t6=lbv;t6<=ubv;t6++) {
              aa[t6][t4] = 0.;;
            }
          }
        }
      }
      if ((t1 >= 1) && (t2 == 0) && (t3 == 0)) {
        for (t4=32*t1;t4<=min(LEN_2D-1,32*t1+31);t4++) {
          dummy(a, b, c, d, e, aa, bb, cc, 0.);;
          for (t5=0;t5<=31;t5++) {
            lbv=0;
            ubv=31;
#pragma ivdep
#pragma vector always
            for (t6=lbv;t6<=ubv;t6++) {
              aa[t6][t4] = 0.;;
            }
          }
        }
      }
      if ((t1 >= 1) && (t2 >= 1) && (t3 == 0)) {
        for (t4=32*t1;t4<=min(LEN_2D-1,32*t1+31);t4++) {
          for (t5=32*t2;t5<=min(100*iterations-1,32*t2+31);t5++) {
            lbv=0;
            ubv=31;
#pragma ivdep
#pragma vector always
            for (t6=lbv;t6<=ubv;t6++) {
              aa[t6][t4] = 0.;;
            }
          }
        }
      }
      if ((t1 >= t3+1) && (t3 >= 1)) {
        for (t4=32*t1;t4<=min(LEN_2D-1,32*t1+31);t4++) {
          for (t5=32*t2;t5<=min(100*iterations-1,32*t2+31);t5++) {
            lbv=32*t3;
            ubv=32*t3+31;
#pragma ivdep
#pragma vector always
            for (t6=lbv;t6<=ubv;t6++) {
              aa[t6][t4] = 0.;;
            }
          }
        }
      }
      if ((t2 == 0) && (t3 == 0)) {
        for (t4=max(LEN_2D,32*t1);t4<=min(100*iterations-1,32*t1+31);t4++) {
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
    
    time_function(&s2102, NULL);
    return EXIT_SUCCESS;
}
