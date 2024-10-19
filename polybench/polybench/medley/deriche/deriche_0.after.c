/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* deriche.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "deriche.h"

#include <omp.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y) ((x) > (y)? (x) : (y))
#define min(x,y) ((x) < (y)? (x) : (y))
#ifndef POLYBENCH_DUMP_ARRAYS
#define DUMP 0
#else
#define DUMP 1
#endif
#ifndef POLYBENCH_CHECKSUM_ARRAYS
#define CHECKSUM 0
#else
#define CHECKSUM 1
#endif


/* Array initialization. */
static
void init_array (int w, int h, DATA_TYPE* alpha,
		 DATA_TYPE POLYBENCH_2D(imgIn,W,H,w,h),
		 DATA_TYPE POLYBENCH_2D(imgOut,W,H,w,h))
{
  int i, j;

  *alpha=0.25; //parameter of the filter

  //input should be between 0 and 1 (grayscale image pixel)
  for (i = 0; i < w; i++)
     for (j = 0; j < h; j++)
	imgIn[i][j] = (DATA_TYPE) ((313*i+991*j)%65536) / 65535.0f;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int w, int h,
		 DATA_TYPE POLYBENCH_2D(imgOut,W,H,w,h))

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("imgOut");
  DATA_TYPE tmp_imgOut = 0;
  for (i = 0; i < w; i++)
    for (j = 0; j < h; j++) {
      if (DUMP) {
        if ((i * h + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
        fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, imgOut[i][j]);
      }
      if (CHECKSUM) tmp_imgOut += imgOut[i][j];
    }
  if (CHECKSUM) {
    fprintf(POLYBENCH_DUMP_TARGET,"\nchecksum: ");
    fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, tmp_imgOut);
  }
  POLYBENCH_DUMP_END("imgOut");
  POLYBENCH_DUMP_FINISH;
}



/* Main computational kernel. The whole function will be timed,
   including the call and return. */
/* Original code provided by Gael Deest */
static
void kernel_deriche(int w, int h, DATA_TYPE alpha,
       DATA_TYPE POLYBENCH_2D(imgIn, W, H, w, h),
       DATA_TYPE POLYBENCH_2D(imgOut, W, H, w, h),
       DATA_TYPE POLYBENCH_2D(y1, W, H, w, h),
       DATA_TYPE POLYBENCH_2D(y2, W, H, w, h)) {
    int i,j;
    DATA_TYPE xm1, tm1, ym1, ym2;
    DATA_TYPE xp1, xp2;
    DATA_TYPE tp1, tp2;
    DATA_TYPE yp1, yp2;

    DATA_TYPE k;
    DATA_TYPE a1, a2, a3, a4, a5, a6, a7, a8;
    DATA_TYPE b1, b2, c1, c2;

#pragma scop
/**/

// Optimized code

// Precompute constants
k = (SCALAR_VAL(1.0) - EXP_FUN(-alpha)) * (SCALAR_VAL(1.0) - EXP_FUN(-alpha)) / (SCALAR_VAL(1.0) + SCALAR_VAL(2.0) * alpha * EXP_FUN(-alpha) - EXP_FUN(SCALAR_VAL(2.0) * alpha));
a1 = k;
a2 = k * EXP_FUN(-alpha) * (alpha - SCALAR_VAL(1.0));
a3 = k * EXP_FUN(-alpha) * (alpha + SCALAR_VAL(1.0));
a4 = -k * EXP_FUN(SCALAR_VAL(-2.0) * alpha);
b1 = POW_FUN(SCALAR_VAL(2.0), -alpha);
b2 = -EXP_FUN(SCALAR_VAL(-2.0) * alpha);
c1 = 1;
c2 = 1;

// Loop unrolling and pipelining for better cache utilization
#pragma omp parallel for private(ym1, ym2, xm1, yp1, yp2, xp1, xp2, tm1, tp1, tp2)
for (int i = 0; i < _PB_W; i++) {
    double ym1 = SCALAR_VAL(0.0);
    double ym2 = SCALAR_VAL(0.0);
    double xm1 = SCALAR_VAL(0.0);
    double yp1 = SCALAR_VAL(0.0);
    double yp2 = SCALAR_VAL(0.0);
    double xp1 = SCALAR_VAL(0.0);
    double xp2 = SCALAR_VAL(0.0);

    // Process the first half of the row
    for (int j = 0; j < _PB_H; j += 2) {
        y1[i][j] = a1 * imgIn[i][j] + a2 * xm1 + b1 * ym1 + b2 * ym2;
        xm1 = imgIn[i][j];
        ym2 = ym1;
        ym1 = y1[i][j];

        if (j + 1 < _PB_H) {
            y1[i][j + 1] = a1 * imgIn[i][j + 1] + a2 * imgIn[i][j] + b1 * y1[i][j] + b2 * ym2;
            xm1 = imgIn[i][j + 1];
            ym2 = y1[i][j];
            ym1 = y1[i][j + 1];
        }
    }

    // Process the second half of the row
    for (int j = _PB_H - 1; j >= 0; j -= 2) {
        y2[i][j] = a3 * xp1 + a4 * xp2 + b1 * yp1 + b2 * yp2;
        xp2 = xp1;
        xp1 = imgIn[i][j];
        yp2 = yp1;
        yp1 = y2[i][j];

        if (j - 1 >= 0) {
            y2[i][j - 1] = a3 * imgIn[i][j] + a4 * xp2 + b1 * y2[i][j] + b2 * yp2;
            xp2 = imgIn[i][j];
            xp1 = imgIn[i][j - 1];
            yp2 = y2[i][j];
            yp1 = y2[i][j - 1];
        }
    }
}

// Combine results
#pragma omp parallel for
for (int i = 0; i < _PB_W; i++) {
    for (int j = 0; j < _PB_H; j++) {
        imgOut[i][j] = c1 * (y1[i][j] + y2[i][j]);
    }
}

// Process columns
#pragma omp parallel for private(tm1, ym1, ym2, tp1, tp2, yp1, yp2)
for (int j = 0; j < _PB_H; j++) {
    double tm1 = SCALAR_VAL(0.0);
    double ym1 = SCALAR_VAL(0.0);
    double ym2 = SCALAR_VAL(0.0);
    double tp1 = SCALAR_VAL(0.0);
    double tp2 = SCALAR_VAL(0.0);
    double yp1 = SCALAR_VAL(0.0);
    double yp2 = SCALAR_VAL(0.0);

    // Process the first half of the column
    for (int i = 0; i < _PB_W; i += 2) {
        y1[i][j] = a5 * imgOut[i][j] + a6 * tm1 + b1 * ym1 + b2 * ym2;
        tm1 = imgOut[i][j];
        ym2 = ym1;
        ym1 = y1[i][j];

        if (i + 1 < _PB_W) {
            y1[i + 1][j] = a5 * imgOut[i + 1][j] + a6 * imgOut[i][j] + b1 * y1[i][j] + b2 * ym2;
            tm1 = imgOut[i + 1][j];
            ym2 = y1[i][j];
            ym1 = y1[i + 1][j];
        }
    }

    // Process the second half of the column
    for (int i = _PB_W - 1; i >= 0; i -= 2) {
        y2[i][j] = a7 * tp1 + a8 * tp2 + b1 * yp1 + b2 * yp2;
        tp2 = tp1;
        tp1 = imgOut[i][j];
        yp2 = yp1;
        yp1 = y2[i][j];

        if (i - 1 >= 0) {
            y2[i - 1][j] = a7 * imgOut[i][j] + a8 * tp2 + b1 * y2[i][j] + b2 * yp2;
            tp2 = imgOut[i][j];
            tp1 = imgOut[i - 1][j];
            yp2 = y2[i][j];
            yp1 = y2[i - 1][j];
        }
    }
}

// Combine results
#pragma omp parallel for
for (int i = 0; i < _PB_W; i++) {
    for (int j = 0; j < _PB_H; j++) {
        imgOut[i][j] = c2 * (y1[i][j] + y2[i][j]);
    }
}
#pragma endscop
}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int w = W;
  int h = H;

  /* Variable declaration/allocation. */
  DATA_TYPE alpha;
  POLYBENCH_2D_ARRAY_DECL(imgIn, DATA_TYPE, W, H, w, h);
  POLYBENCH_2D_ARRAY_DECL(imgOut, DATA_TYPE, W, H, w, h);
  POLYBENCH_2D_ARRAY_DECL(y1, DATA_TYPE, W, H, w, h);
  POLYBENCH_2D_ARRAY_DECL(y2, DATA_TYPE, W, H, w, h);


  /* Initialize array(s). */
  init_array (w, h, &alpha, POLYBENCH_ARRAY(imgIn), POLYBENCH_ARRAY(imgOut));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_deriche (w, h, alpha, POLYBENCH_ARRAY(imgIn), POLYBENCH_ARRAY(imgOut), POLYBENCH_ARRAY(y1), POLYBENCH_ARRAY(y2));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(w, h, POLYBENCH_ARRAY(imgOut)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(imgIn);
  POLYBENCH_FREE_ARRAY(imgOut);
  POLYBENCH_FREE_ARRAY(y1);
  POLYBENCH_FREE_ARRAY(y2);

  return 0;
}