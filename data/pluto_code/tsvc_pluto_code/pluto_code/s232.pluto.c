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


real_t s232(struct args_t * func_args)
{

//    loop interchange
//    interchanging of triangular loops

    initialise_arrays(__func__);
    func_args->c1 = omp_get_wtime();
  int t1, t2, t3, t4, t5, t6;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
lbp=0;
ubp=floord(100*iterations+LEN_2D-260,32);
#pragma omp parallel for private(lbv,ubv,t2,t3,t4,t5,t6)
for (t1=lbp;t1<=ubp;t1++) {
  for (t2=0;t2<=min(min(min(min(floord(25*iterations-1,8),floord(-2500024*t1+7812375*iterations+19687697,2499960)),floord(-19999992*t1+62499975*iterations-2499999,19999480)),floord(-26666656*t1+83333300*iterations+833333*LEN_2D-216666580,26665976)),floord(-80000768*t1+249996025*iterations+2500024*LEN_2D-12499865,79998728));t2++) {
    for (t3=t2;t3<=min(min(min(min(min(min(min(min(floord(-106669376*t1+333341800*iterations+6666665*LEN_2D-1719999912,106663904),floord(-320008160*t1+1000025500*iterations+10000255*LEN_2D-2600066300,319991712)),floord(-9184*t1+319991712*t2+28700*iterations+287*LEN_2D-74620,319991712)),floord(-9152*t1+319991712*t2+28600*iterations+10000027*LEN_2D-2560008056,319991712)),floord(-320011360*t1+999984100*iterations+10000355*LEN_2D+2539959100,319994912)),floord(-320011328*t1+999984100*iterations+20000195*LEN_2D-30000036,319994912)),floord(-9184*t1+319994912*t2+287*LEN_2D+2869954080,319994912)),floord(-9152*t1+319994912*t2+10000127*LEN_2D+299994944,319994912)),t1+t2+1);t3++) {
      if (t3 >= 1) {
        for (t4=max(max(max(1,32*t1),-32*t2+32*t3-31),32*t3-100*iterations+1);t4<=min(LEN_2D-1,32*t1+31);t4++) {
          for (t5=max(32*t2,32*t3-t4);t5<=min(min(100*iterations-1,32*t2+31),32*t3+30);t5++) {
            for (t6=max(32*t3,t5+1);t6<=min(32*t3+31,t4+t5);t6++) {
              aa[t4][(-t5+t6)] = aa[t4][(-t5+t6)-1]*aa[t4][(-t5+t6)-1]+bb[t4][(-t5+t6)];;
            }
          }
        }
      }
      if ((t1 == 0) && (t2 == 0) && (t3 == 0)) {
        dummy(a, b, c, d, e, aa, bb, cc, 1.);;
      }
      if ((t2 == 0) && (t3 == 0)) {
        for (t4=max(1,32*t1);t4<=min(LEN_2D-1,32*t1+31);t4++) {
          dummy(a, b, c, d, e, aa, bb, cc, 1.);;
          for (t5=0;t5<=30;t5++) {
            for (t6=t5+1;t6<=min(31,t4+t5);t6++) {
              aa[t4][(-t5+t6)] = aa[t4][(-t5+t6)-1]*aa[t4][(-t5+t6)-1]+bb[t4][(-t5+t6)];;
            }
          }
        }
      }
      if ((t2 == 0) && (t3 == 0)) {
        for (t4=max(LEN_2D,32*t1);t4<=min(100*iterations-1,32*t1+31);t4++) {
          dummy(a, b, c, d, e, aa, bb, cc, 1.);;
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
    
    time_function(&s232, NULL);
    return EXIT_SUCCESS;
}
