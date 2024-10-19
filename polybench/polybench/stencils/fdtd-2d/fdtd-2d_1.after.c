/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* fdtd-2d.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "fdtd-2d.h"

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
void init_array (int tmax,
		 int nx,
		 int ny,
		 DATA_TYPE POLYBENCH_2D(ex,NX,NY,nx,ny),
		 DATA_TYPE POLYBENCH_2D(ey,NX,NY,nx,ny),
		 DATA_TYPE POLYBENCH_2D(hz,NX,NY,nx,ny),
		 DATA_TYPE POLYBENCH_1D(_fict_,TMAX,tmax))
{
  int i, j;

  for (i = 0; i < tmax; i++)
    _fict_[i] = (DATA_TYPE) i;
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      {
	ex[i][j] = ((DATA_TYPE) i*(j+1)) / nx;
	ey[i][j] = ((DATA_TYPE) i*(j+2)) / ny;
	hz[i][j] = ((DATA_TYPE) i*(j+3)) / nx;
      }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int nx,
		 int ny,
		 DATA_TYPE POLYBENCH_2D(ex,NX,NY,nx,ny),
		 DATA_TYPE POLYBENCH_2D(ey,NX,NY,nx,ny),
		 DATA_TYPE POLYBENCH_2D(hz,NX,NY,nx,ny))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("ex");
  DATA_TYPE tmp_ex = 0;
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++) {
      if (DUMP) {
        if ((i * nx + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
        fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, ex[i][j]);
      }
      if (CHECKSUM) tmp_ex += ex[i][j];
    }
  if (CHECKSUM) {
    fprintf(POLYBENCH_DUMP_TARGET,"\nchecksum: ");
    fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, tmp_ex);
  }
  POLYBENCH_DUMP_END("ex");
  POLYBENCH_DUMP_FINISH;

  POLYBENCH_DUMP_BEGIN("ey");
  DATA_TYPE tmp_ey = 0;
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++) {
      if (DUMP) {
        if ((i * nx + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
        fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, ey[i][j]);
      }
      if (CHECKSUM) tmp_ey += ey[i][j];
    }
  if (CHECKSUM) {
    fprintf(POLYBENCH_DUMP_TARGET,"\nchecksum: ");
    fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, tmp_ey);
  }
  POLYBENCH_DUMP_END("ey");

  POLYBENCH_DUMP_BEGIN("hz");
  DATA_TYPE tmp_hz = 0;
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++) {
      if (DUMP) {
        if ((i * nx + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
        fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, hz[i][j]);
      }
      if (CHECKSUM) tmp_hz += hz[i][j];
    }
  if (CHECKSUM) {
    fprintf(POLYBENCH_DUMP_TARGET,"\nchecksum: ");
    fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, tmp_hz);
  }
  POLYBENCH_DUMP_END("hz");
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_fdtd_2d(int tmax,
		    int nx,
		    int ny,
		    DATA_TYPE POLYBENCH_2D(ex,NX,NY,nx,ny),
		    DATA_TYPE POLYBENCH_2D(ey,NX,NY,nx,ny),
		    DATA_TYPE POLYBENCH_2D(hz,NX,NY,nx,ny),
		    DATA_TYPE POLYBENCH_1D(_fict_,TMAX,tmax))
{
  int t, i, j;

#pragma scop
/**/

int t1, t2, t3, t4, t5, t6, t7, t8;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

lbp = 0;
ubp = floord(_PB_TMAX - 1, 32);
#pragma omp parallel for private(lbv, ubv, t3, t4, t5, t6, t7, t8)
for (t2 = lbp; t2 <= ubp; t2++) {
    for (t3 = 0; t3 <= floord(_PB_NY - 1, 32); t3++) {
        for (t4 = max(0, 32 * t2); t4 <= min(_PB_TMAX - 1, 32 * t2 + 31); t4++) {
            for (t5 = max(0, 32 * t3); t5 <= min(_PB_NY - 1, 32 * t3 + 31); t5++) {
                ey[0][t5] = _fict_[t4];
            }
        }
    }
}

lbp = 0;
ubp = floord(_PB_NX - 1, 32);
#pragma omp parallel for private(lbv, ubv, t3, t4, t5, t6, t7, t8)
for (t2 = lbp; t2 <= ubp; t2++) {
    for (t3 = 0; t3 <= floord(_PB_NY - 1, 32); t3++) {
        for (t4 = max(1, 32 * t2); t4 <= min(_PB_NX - 1, 32 * t2 + 31); t4++) {
            for (t5 = max(0, 32 * t3); t5 <= min(_PB_NY - 1, 32 * t3 + 31); t5++) {
                ey[t4][t5] = ey[t4][t5] - SCALAR_VAL(0.5) * (hz[t4][t5] - hz[t4 - 1][t5]);
            }
        }
    }
}

lbp = 0;
ubp = floord(_PB_NX - 1, 32);
#pragma omp parallel for private(lbv, ubv, t3, t4, t5, t6, t7, t8)
for (t2 = lbp; t2 <= ubp; t2++) {
    for (t3 = 0; t3 <= floord(_PB_NY - 1, 32); t3++) {
        for (t4 = max(0, 32 * t2); t4 <= min(_PB_NX - 1, 32 * t2 + 31); t4++) {
            for (t5 = max(1, 32 * t3); t5 <= min(_PB_NY - 1, 32 * t3 + 31); t5++) {
                ex[t4][t5] = ex[t4][t5] - SCALAR_VAL(0.5) * (hz[t4][t5] - hz[t4][t5 - 1]);
            }
        }
    }
}

lbp = 0;
ubp = floord(_PB_NX - 2, 32);
#pragma omp parallel for private(lbv, ubv, t3, t4, t5, t6, t7, t8)
for (t2 = lbp; t2 <= ubp; t2++) {
    for (t3 = 0; t3 <= floord(_PB_NY - 2, 32); t3++) {
        for (t4 = max(0, 32 * t2); t4 <= min(_PB_NX - 2, 32 * t2 + 31); t4++) {
            for (t5 = max(0, 32 * t3); t5 <= min(_PB_NY - 2, 32 * t3 + 31); t5++) {
                hz[t4][t5] = hz[t4][t5] - SCALAR_VAL(0.7) * (ex[t4][t5 + 1] - ex[t4][t5] + ey[t4 + 1][t5] - ey[t4][t5]);
            }
        }
    }
}
#pragma endscop
}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int tmax = TMAX;
  int nx = NX;
  int ny = NY;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(ex,DATA_TYPE,NX,NY,nx,ny);
  POLYBENCH_2D_ARRAY_DECL(ey,DATA_TYPE,NX,NY,nx,ny);
  POLYBENCH_2D_ARRAY_DECL(hz,DATA_TYPE,NX,NY,nx,ny);
  POLYBENCH_1D_ARRAY_DECL(_fict_,DATA_TYPE,TMAX,tmax);

  /* Initialize array(s). */
  init_array (tmax, nx, ny,
	      POLYBENCH_ARRAY(ex),
	      POLYBENCH_ARRAY(ey),
	      POLYBENCH_ARRAY(hz),
	      POLYBENCH_ARRAY(_fict_));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_fdtd_2d (tmax, nx, ny,
		  POLYBENCH_ARRAY(ex),
		  POLYBENCH_ARRAY(ey),
		  POLYBENCH_ARRAY(hz),
		  POLYBENCH_ARRAY(_fict_));


  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(nx, ny, POLYBENCH_ARRAY(ex),
				    POLYBENCH_ARRAY(ey),
				    POLYBENCH_ARRAY(hz)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(ex);
  POLYBENCH_FREE_ARRAY(ey);
  POLYBENCH_FREE_ARRAY(hz);
  POLYBENCH_FREE_ARRAY(_fict_);

  return 0;
}
