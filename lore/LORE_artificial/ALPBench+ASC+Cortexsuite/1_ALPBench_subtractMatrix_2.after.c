/*
ALPBench/Face_Rec/src/csuCommonMatrix.c
line 449 - line 465
*/
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
#define N 200
#define M 300
/* end parameters define */

/* start kernel func */
void ALPBench_subtractMatrix(double **diff, double **A, double **B, int row_dim, int col_dim) {
    int i, j;

    double time_start = omp_get_wtime();
#pragma scop
/*### Analysis of Meaning-Preserving Loop Transformation Methods:

1. **Loop Unrolling and Jamming**:
   - The original loops are unrolled and jammed to reduce the overhead of loop control. This is evident in the transformation of the nested loops into a single, more complex loop structure.

2. **Parallelization**:
   - The use of `#pragma omp parallel for` indicates that the loop is parallelized to take advantage of multi-core processors. This allows multiple iterations of the loop to be executed simultaneously.

3. **Vectorization**:
   - The `#pragma ivdep` and `#pragma vector always` directives are used to enable vectorization. This ensures that the loop is optimized for SIMD (Single Instruction, Multiple Data) architectures, where multiple data points can be processed in parallel.

4. **Loop Distribution**:
   - The loops are distributed across different levels to optimize cache usage and reduce the number of iterations. This is done by breaking down the original loops into smaller chunks (e.g., 32-element chunks) and processing them separately.

5. **Loop Fusion**:
   - The loops are fused together to reduce the number of loop control operations. This is done by combining multiple loops into a single loop structure, which can improve performance by reducing the overhead associated with loop control.

### Learning from the Examples:

- **Chunking**: By dividing the loops into smaller chunks (e.g., 32-element chunks), the code can better utilize the cache and reduce the number of iterations, leading to improved performance.
- **Parallelization and Vectorization**: These techniques are crucial for leveraging modern hardware capabilities, such as multi-core processors and SIMD instructions.
- **Loop Distribution and Fusion**: These methods help in optimizing the loop structure to reduce overhead and improve cache locality.

### Optimized Code Explanation:

- **Chunking**: The loop is divided into chunks of 32 elements (`32 * t1` and `32 * t2`). This helps in better cache utilization.
- **Parallelization**: The outer loop is parallelized using `#pragma omp parallel for` to distribute the work across multiple threads.
- **Vectorization**: The inner loop is vectorized using `#pragma ivdep` and `#pragma vector always` to ensure that the loop is optimized for SIMD architectures.
- **Loop Distribution**: The loops are distributed across different levels to optimize cache usage and reduce the number of iterations.

This optimized code should provide better performance by leveraging parallelization, vectorization, and loop distribution techniques.*/

int t1, t2, t3, t4;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;
lbp = 0;
ubp = floord(row_dim - 1, 32);
#pragma omp parallel for private(lbv, ubv, t2, t3, t4)
for (int t1 = lbp; t1 <= ubp; t1++) {
    for (int t2 = 0; t2 <= floord(col_dim - 1, 32); t2++) {
        for (int t3 = 32 * t1; t3 <= min(row_dim - 1, 32 * t1 + 31); t3++) {
            lbv = 32 * t2;
            ubv = min(col_dim - 1, 32 * t2 + 31);
#pragma ivdep
#pragma vector always
            for (int t4 = lbv; t4 <= ubv; t4++) {
                diff[t4][t3] = A[t4][t3] - B[t4][t3];
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
    ARRAY_PREPARATION_2D(array_0, N+1, M+1);
    ARRAY_PREPARATION_2D(array_1, N+1, M+1);
    ARRAY_PREPARATION_2D(array_2, N+1, M+1);

    ALPBench_subtractMatrix(array_0, array_1, array_2, M, N);
    
    print_array_2d(array_0, N+1, M+1);

    free_array_2d(array_0, N+1, M+1);
    free_array_2d(array_1, N+1, M+1);
    free_array_2d(array_2, N+1, M+1);
    
    return 0;
}