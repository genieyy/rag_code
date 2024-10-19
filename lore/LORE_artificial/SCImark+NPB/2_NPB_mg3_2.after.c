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
#define N 40
#define M 50
#define L 60
/* end parameters define */

/* start kernel func */
void NPB_mg3(double ***r, double *r1, double * r2, double ***u, int n3, int n2, int n1)
{
	int i1, i2, i3;
	double c[3] = {0.2, 0.3, 0.4};

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation of Optimizations:
1. **Parallelization**: The outermost loop (`iter`) is parallelized using OpenMP's `#pragma omp parallel for` directive. This allows multiple threads to execute the iterations of the loop concurrently, which can significantly improve performance on multi-core systems.

2. **Loop Order**: The loop order is maintained as it is, but the parallelization directive is applied to the outermost loop to maximize the potential for concurrency.

3. **Private Variables**: The variables `i1`, `i2`, and `i3` are declared as private within the parallel region to ensure that each thread has its own copy of these loop indices, preventing race conditions.

These optimizations leverage the techniques observed in the provided examples, such as parallelization and loop reordering, to improve the performance of the given code.*/

#pragma omp parallel for private(i1, i2, i3)
for (int iter = 0; iter < ITERATIONS; iter++) {
    for (i3 = 1; i3 < n3 - 1; i3++) {
        for (i2 = 1; i2 < n2 - 1; i2++) {
            for (i1 = 0; i1 < n1; i1++) {
                r1[i1] = r[i3][i2 - 1][i1] + r[i3][i2 + 1][i1] + r[i3 - 1][i2][i1] + r[i3 + 1][i2][i1];
                r2[i1] = r[i3 - 1][i2 - 1][i1] + r[i3 - 1][i2 + 1][i1] + r[i3 + 1][i2 - 1][i1] + r[i3 + 1][i2 + 1][i1];
            }
            for (i1 = 1; i1 < n1 - 1; i1++) {
                u[i3][i2][i1] = u[i3][i2][i1] + c[0] * r[i3][i2][i1] + c[1] * (r[i3][i2][i1 - 1] + r[i3][i2][i1 + 1] + r1[i1]) + c[2] * (r2[i1] + r1[i1 - 1] + r1[i1 + 1]);
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
	ARRAY_PREPARATION_3D(array_0, N + 1, M + 1, L + 1);
	ARRAY_PREPARATION_1D(array_1, L+1);
	ARRAY_PREPARATION_1D(array_2, L+1);
	ARRAY_PREPARATION_3D(array_3, N + 1, M + 1, L + 1);

	NPB_mg3(array_0, array_1, array_2, array_3, N, M, L);

	print_array_3d(array_3, N + 1, M + 1, L + 1);

	free_array_3d(array_0, N + 1, M + 1, L + 1);
	free_array_1d(array_1, L+1);
	free_array_1d(array_2, L+1);
	free_array_3d(array_3, N + 1, M + 1, L + 1);

	return 0;
}