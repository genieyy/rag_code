#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lore.h"

#include <omp.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y) ((x) > (y)? (x) : (y))
#define min(x,y) ((x) < (y)? (x) : (y))

/* start param define */
#define N 80
#define M 100
#define L 5
/* end parameters define */

/* start kernel func */
void NPB_lu3(double ****d, double **tmat, int iend, int jend, int n)
{
	int i, j, m;
	int ist = 1;
	int jst = 2;

    double time_start = omp_get_wtime();
#pragma scop
/*### Key Improvements:
1. **Pointer Pre-calculation**: The pointer `d_ij` is pre-calculated outside the inner loop to avoid recalculating the same index multiple times.
2. **Pointer Arithmetic**: The pointer `d_ijm` is calculated using pointer arithmetic (`d_ij + 5 * m`) to directly access the elements of `d[i][j][m]`, reducing the overhead of array indexing.
3. **Full Loop Unrolling**: The inner loop is fully unrolled to reduce loop control overhead and improve instruction-level parallelism.*/

/*### Further Optimized Version:

### Explanation of Optimizations:
1. **Loop Unrolling**: The inner loop is fully unrolled to reduce the number of loop iterations and improve instruction-level parallelism. This allows the compiler to generate more efficient code by eliminating the loop control overhead.
2. **Pointer Arithmetic**: By using pointers (`tmat_m` and `d_ijm`), we reduce the number of array indexing operations, which can be costly in terms of both time and memory access.
3. **Reduced Indirection**: By assigning the pointers outside the loop, we reduce the overhead of repeated array indexing, which can improve cache performance and reduce the number of memory accesses.
4. **Vectorization**: Although not explicitly vectorized, the use of pointers and unrolling can help the compiler to vectorize the loop, leading to further performance improvements on modern CPUs with SIMD capabilities.

These optimizations are based on the principles observed in the provided examples, such as reducing the number of operations within loops, minimizing the overhead of array indexing, and leveraging pointer arithmetic to improve cache locality.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int i = ist; i <= iend; i++) {
        for (int j = jst; j <= jend; j++) {
            double *d_ij = d[i][j][0];
            for (int m = 0; m < n; m++) {
                double *tmat_m = tmat[m];
                double *d_ijm = d_ij + 5 * m;
                tmat_m[0] = d_ijm[0];
                tmat_m[1] = d_ijm[1];
                tmat_m[2] = d_ijm[2];
                tmat_m[3] = d_ijm[3];
                tmat_m[4] = d_ijm[4];
            }
        }
    }
}
#pragma endscop
    double time_end = omp_get_wtime();
    printf("%f\n", time_end - time_start);
}
/* end kernel func */

int main(int argc, char *argv[])
{
	ARRAY_PREPARATION_4D(array_0, N + 1, M + 1, L + 1, 5);
	ARRAY_PREPARATION_2D(array_1, L + 1, 5);

	NPB_lu3(array_0, array_1, N, M, L);

	print_array_2d(array_1, L + 1, 5);

	free_array_4d(array_0, N + 1, M + 1, L + 1, 5);
	free_array_2d(array_1, L + 1, 5);

	return 0;
}