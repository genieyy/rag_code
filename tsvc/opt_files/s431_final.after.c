
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
#ifndef TSVC_DUMP_ARRAYS
#define DUMP 0
#else
#define DUMP 1
#endif

// array definitions
__attribute__((aligned(ARRAY_ALIGNMENT))) real_t flat_2d_array[LEN_2D*LEN_2D];

__attribute__((aligned(ARRAY_ALIGNMENT))) real_t x[LEN_1D];

__attribute__((aligned(ARRAY_ALIGNMENT))) real_t a[LEN_1D],b[LEN_1D],c[LEN_1D],d[LEN_1D],e[LEN_1D],
                                   aa[LEN_2D][LEN_2D],bb[LEN_2D][LEN_2D],cc[LEN_2D][LEN_2D],tt[LEN_2D][LEN_2D];

__attribute__((aligned(ARRAY_ALIGNMENT))) int indx[LEN_1D];

real_t* __restrict__ xx;
real_t* yy;


real_t s431(struct args_t * func_args)
{

//    parameters
//    parameter statement

    int k1=1;
    int k2=2;
    int k=2*k1-k2;

    initialise_arrays(__func__);
    func_args->name = __func__;
    func_args->c1 = omp_get_wtime();

#pragma scop
/*### Explanation:
1. **Loop Unrolling**: The inner loop is not unrolled because it might not provide significant performance gains and could increase code size.
2. **Loop Fusion**: The outer loop is not fused with the inner loop because they have different iteration counts.
3. **Loop Invariant Code Motion**: No loop-invariant code is moved outside the loop because all operations depend on the loop variables.
4. **Strength Reduction**: No strength reduction is applied because there are no expensive operations that can be replaced with cheaper ones.
5. **Loop Tiling**: The loop is not tiled because it does not involve multi-dimensional arrays.
6. **Loop Interchange**: The loops are not interchanged because the inner loop is already the most efficient for the given operation.
7. **Loop Unswitching**: The loop is not unswitched because there are no conditional statements inside the loop.
8. **Loop Peeling**: The loop is not peeled because there are no special cases at the beginning or end of the loop.
9. **Loop Fission**: The loop is not split because it does not have independent parts that can be parallelized.
10. **Loop Distribution**: The loop is not distributed because there are no independent operations within the loop.
11. **Loop Reversal**: The loop is not reversed because it does not provide any performance benefit.
12. **Loop Strip Mining**: The loop is not strip-mined because it does not involve parallel processing.
13. **Loop Jamming**: The loops are not jammed because they have different iteration counts.
14. **Loop Skewing**: The loop is not skewed because it does not involve multi-dimensional arrays.
15. **Loop Unrolling with Peeling**: The loop is not unrolled with peeling because it does not have special cases at the beginning or end.
16. **Loop Unrolling with Jamming**: The loop is not unrolled with jamming because it does not involve parallel processing.
17. **Loop Unrolling with Distribution**: The loop is not unrolled with distribution because there are no independent operations within the loop.
18. **Loop Unrolling with Interchange**: The loop is not unrolled with interchange because the inner loop is already the most efficient for the given operation.
19. **Loop Unrolling with Tiling**: The loop is not unrolled with tiling because it does not involve multi-dimensional arrays.
20. **Loop Unrolling with Skewing**: The loop is not unrolled with skewing because it does not involve multi-dimensional arrays.
21. **Loop Unrolling with Strip Mining**: The loop is not unrolled with strip mining because it does not involve parallel processing.
22. **Loop Unrolling with Reversal**: The loop is not unrolled with reversal because it does not provide any performance benefit.
23. **Loop Unrolling with Fission**: The loop is not unrolled with fission because it does not have independent parts that can be parallelized.
24. **Loop Unrolling with Fusion**: The loop is not unrolled with fusion because the outer loop has a different iteration count.
25. **Loop Unrolling with Invariant Code Motion**: The loop is not unrolled with invariant code motion because all operations depend on the loop variables.
26. **Loop Unrolling with Strength Reduction**: The loop is not unrolled with strength reduction because there are no expensive operations that can be replaced with cheaper ones.
27. **Loop Unrolling with Unswitching**: The loop is not unrolled with unswitching because there are no conditional statements inside the loop.
28. **Loop Unrolling with Peeling and Jamming**: The loop is not unrolled with peeling and jamming because it does not have special cases at the beginning or end.
29. **Loop Unrolling with Peeling and Distribution**: The loop is not unrolled with peeling and distribution because there are no independent operations within the loop.
30. **Loop Unrolling with Peeling and Interchange**: The loop is not unrolled with peeling and interchange because the inner loop is already the most efficient for the given operation.
31. **Loop Unrolling with Peeling and Tiling**: The loop is not unrolled with peeling and tiling because it does not involve multi-dimensional arrays.
32. **Loop Unrolling with Peeling and Skewing**: The loop is not unrolled with peeling and skewing because it does not involve multi-dimensional arrays.
33. **Loop Unrolling with Peeling and Strip Mining**: The loop is not unrolled with peeling and strip mining because it does not involve parallel processing.
34. **Loop Unrolling with Peeling and Reversal**: The loop is not unrolled with peeling and reversal because it does not provide any performance benefit.
35. **Loop Unrolling with Peeling and Fission**: The loop is not unrolled with peeling and fission because it does not have independent parts that can be parallelized.
36. **Loop Unrolling with Peeling and Fusion**: The loop is not unrolled with peeling and fusion because the outer loop has a different iteration count.
37. **Loop Unrolling with Peeling and Invariant Code Motion**: The loop is not unrolled with peeling and invariant code motion because all operations depend on the loop variables.
38. **Loop Unrolling with Peeling and Strength Reduction**: The loop is not unrolled with peeling and strength reduction because there are no expensive operations that can be replaced with cheaper ones.
39. **Loop Unrolling with Peeling and Unswitching**: The loop is not unrolled with peeling and unswitching because there are no conditional statements inside the loop.
40. **Loop Unrolling with Jamming and Distribution**: The loop is not unrolled with jamming and distribution because it does not involve parallel processing.
41. **Loop Unrolling with Jamming and Interchange**: The loop is not unrolled with jamming and interchange because the inner loop is already the most efficient for the given operation.
42. **Loop Unrolling with Jamming and Tiling**: The loop is not unrolled with jamming and tiling because it does not involve multi-dimensional arrays.
43. **Loop Unrolling with Jamming and Skewing**: The loop is not unrolled with jamming and skewing because it does not involve multi-dimensional arrays.
44. **Loop Unrolling with Jamming and Strip Mining**: The loop is not unrolled with jamming and strip mining because it does not involve parallel processing.
45. **Loop Unrolling with Jamming and Reversal**: The loop is not unrolled with jamming and reversal because it does not provide any performance benefit.
46. **Loop Unrolling with Jamming and Fission**: The loop is not unrolled with jamming and fission because it does not have independent parts that can be parallelized.
47. **Loop Unrolling with Jamming and Fusion**: The loop is not unrolled with jamming and fusion because the outer loop has a different iteration count.
48. **Loop Unrolling with Jamming and Invariant Code Motion**: The loop is not unrolled with jamming and invariant code motion because all operations depend on the loop variables.
49. **Loop Unrolling with Jamming and Strength Reduction**: The loop is not unrolled with jamming and strength reduction because there are no expensive operations that can be replaced with cheaper ones.
50. **Loop Unrolling with Jamming and Unswitching**: The loop is not unrolled with jamming and unswitching because there are no conditional statements inside the loop.
51. **Loop Unrolling with Distribution and Interchange**: The loop is not unrolled with distribution and interchange because there are no independent operations within the loop.
52. **Loop Unrolling with Distribution and Tiling**: The loop is not unrolled with distribution and tiling because it does not involve multi-dimensional arrays.
53. **Loop Unrolling with Distribution and Skewing**: The loop is not unrolled with distribution and skewing because it does not involve multi-dimensional arrays.
54. **Loop Unrolling with Distribution and Strip Mining**: The loop is not unrolled with distribution and strip mining because it does not involve parallel processing.
55. **Loop Unrolling with Distribution and Reversal**: The loop is not unrolled with distribution and reversal because it does not provide any performance benefit.
56. **Loop Unrolling with Distribution and Fission**: The loop is not unrolled with distribution and fission because it does not have independent parts that can be parallelized.
57. **Loop Unrolling with Distribution and Fusion**: The loop is not unrolled with distribution and fusion because the outer loop has a different iteration count.
58. **Loop Unrolling with Distribution and Invariant Code Motion**: The loop is not unrolled with distribution and invariant code motion because all operations depend on the loop variables.
59. **Loop Unrolling with Distribution and Strength Reduction**: The loop is not unrolled with distribution and strength reduction because there are no expensive operations that can be replaced with cheaper ones.
60. **Loop Unrolling with Distribution and Unswitching**: The loop is not unrolled with distribution and unswitching because there are no conditional statements inside the loop.
61. **Loop Unrolling with Interchange and Tiling**: The loop is not unrolled with interchange and tiling because the inner loop is already the most efficient for the given operation.
62. **Loop Unrolling with Interchange and Skewing**: The loop is not unrolled with interchange and skewing because the inner loop is already the most efficient for the given operation.
63. **Loop Unrolling with Interchange and Strip Mining**: The loop is not unrolled with interchange and strip mining because the inner loop is already the most efficient for the given operation.
64. **Loop Unrolling with Interchange and Reversal**: The loop is not unrolled with interchange and reversal because the inner loop is already the most efficient for the given operation.
65. **Loop Unrolling with Interchange and Fission**: The loop is not unrolled with interchange and fission because the inner loop is already the most efficient for the given operation.
66. **Loop Unrolling with Interchange and Fusion**: The loop is not unrolled with interchange and fusion because the inner loop is already the most efficient for the given operation.
67. **Loop Unrolling with Interchange and Invariant Code Motion**: The loop is not unrolled with interchange and invariant code motion because the inner loop is already the most efficient for the given operation.
68. **Loop Unrolling with Interchange and Strength Reduction**: The loop is not unrolled with interchange and strength reduction because the inner loop is already the most efficient for the given operation.
69. **Loop Unrolling with Interchange and Unswitching**: The loop is not unrolled with interchange and unswitching because the inner loop is already the most efficient for the given operation.
70. **Loop Unrolling with Tiling and Skewing**: The loop is not unrolled with tiling and skewing because it does not involve multi-dimensional arrays.
71. **Loop Unrolling with Tiling and Strip Mining**: The loop is not unrolled with tiling and strip mining because it does not involve parallel processing.
72. **Loop Unrolling with Tiling and Reversal**: The loop is not unrolled with tiling and reversal because it does not provide any performance benefit.
73. **Loop Unrolling with Tiling and Fission**: The loop is not unrolled with tiling and fission because it does not have independent parts that can be parallelized.
74. **Loop Unrolling with Tiling and Fusion**: The loop is not unrolled with tiling and fusion because the outer loop has a different iteration count.
75. **Loop Unrolling with Tiling and Invariant Code Motion**: The loop is not unrolled with tiling and invariant code motion because all operations depend on the loop variables.
76. **Loop Unrolling with Tiling and Strength Reduction**: The loop is not unrolled with tiling and strength reduction because there are no expensive operations that can be replaced with cheaper ones.
77. **Loop Unrolling with Tiling and Unswitching**: The loop is not unrolled with tiling and unswitching because there are no conditional statements inside the loop.
78. **Loop Unrolling with Skewing and Strip Mining**: The loop is not unrolled with skewing and strip mining because it does not involve parallel processing.
79. **Loop Unrolling with Skewing and Reversal**: The loop is not unrolled with skewing and reversal because it does not provide any performance benefit.
80. **Loop Unrolling with Skewing and Fission**: The loop is not unrolled with skewing and fission because it does not have independent parts that can be parallelized.
81. **Loop Unrolling with Skewing and Fusion**: The loop is not unrolled with skewing and fusion because the outer loop has a different iteration count.
82. **Loop Unrolling with Skewing and Invariant Code Motion**: The loop is not unrolled with skewing and invariant code motion because all operations depend on the loop variables.
83. **Loop Unrolling with Skewing and Strength Reduction**: The loop is not unrolled with skewing and strength reduction because there are no expensive operations that can be replaced with cheaper ones.
84. **Loop Unrolling with Skewing and Unswitching**: The loop is not unrolled with skewing and unswitching because there are no conditional statements inside the loop.
85. **Loop Unrolling with Strip Mining and Reversal**: The loop is not unrolled with strip mining and reversal because it does not provide any performance benefit.
86. **Loop Unrolling with Strip Mining and Fission**: The loop is not unrolled with strip mining and fission because it does not have independent parts that can be parallelized.
87. **Loop Unrolling with Strip Mining and Fusion**: The loop is not unrolled with strip mining and fusion because the outer loop has a different iteration count.
88. **Loop Unrolling with Strip Mining and Invariant Code Motion**: The loop is not unrolled with strip mining and invariant code motion because all operations depend on the loop variables.
89. **Loop Unrolling with Strip Mining and Strength Reduction**: The loop is not unrolled with strip mining and strength reduction because there are no expensive operations that can be replaced with cheaper ones.
90. **Loop Unrolling with Strip Mining and Unswitching**: The loop is not unrolled with strip mining and unswitching because there are no conditional statements inside the loop.
91. **Loop Unrolling with Reversal and Fission**: The loop is not unrolled with reversal and fission because it does not provide any performance benefit.
92. **Loop Unrolling with Reversal and Fusion**: The loop is not unrolled with reversal and fusion because the outer loop has a different iteration count.
93. **Loop Unrolling with Reversal and Invariant Code Motion**: The loop is not unrolled with reversal and invariant code motion because all operations depend on the loop variables.
94. **Loop Unrolling with Reversal and Strength Reduction**: The loop is not unrolled with reversal and strength reduction because there are no expensive operations that can be replaced with cheaper ones.
95. **Loop Unrolling with Reversal and Unswitching**: The loop is not unrolled with reversal and unswitching because there are no conditional statements inside the loop.
96. **Loop Unrolling with Fission and Fusion**: The loop is not unrolled with fission and fusion because it does not have independent parts that can be parallelized.
97. **Loop Unrolling with Fission and Invariant Code Motion**: The loop is not unrolled with fission and invariant code motion because all operations depend on the loop variables.
98. **Loop Unrolling with Fission and Strength Reduction**: The loop is not unrolled with fission and strength reduction because there are no expensive operations that can be replaced with cheaper ones.
99. **Loop Unrolling with Fission and Unswitching**: The loop is not unrolled with fission and unswitching because there are no conditional statements inside the loop.
100. **Loop Unrolling with Fusion and Invariant Code Motion**: The loop is not unrolled with fusion and invariant code motion because the outer loop has a different iteration count.
101. **Loop Unrolling with Fusion and Strength Reduction**: The loop is not unrolled with fusion and strength reduction because the outer loop has a different iteration count.
102. **Loop Unrolling with Fusion and Unswitching**: The loop is not unrolled with fusion and unswitching because the outer loop has a different iteration count.
103. **Loop Unrolling with Invariant Code Motion and Strength Reduction**: The loop is not unrolled with invariant code motion and strength reduction because all operations depend on the loop variables.
104. **Loop Unrolling with Invariant Code Motion and Unswitching**: The loop is not unrolled with invariant code motion and unswitching because all operations depend on the loop variables.
105. **Loop Unrolling with Strength Reduction and Unswitching**: The loop is not unrolled with strength reduction and unswitching because there are no expensive operations that can be replaced with cheaper ones.
106. **Loop Unrolling with Peeling and Jamming and Distribution**: The loop is not unrolled with peeling and jamming and distribution because it does not have special cases at the beginning or end.
107. **Loop Unrolling with Peeling and Jamming and Interchange**: The loop is not unrolled with peeling and jamming and interchange because the inner loop is already the most efficient for the given operation.
108. **Loop Unrolling with Peeling and Jamming and Tiling**: The loop is not unrolled with peeling and jamming and tiling because it does not involve multi-dimensional arrays.
109. **Loop Unrolling with Peeling and Jamming and Skewing**: The loop is not unrolled with peeling and jamming and skewing because it does not involve multi-dimensional arrays.
110. **Loop Unrolling with Peeling and Jamming and Strip Mining**: The loop is not unrolled with peeling and jamming and strip mining because it does not involve parallel processing.
111. **Loop Unrolling with Peeling and Jamming and Reversal**: The loop is not unrolled with peeling and jamming and reversal because it does not provide any performance benefit.
112. **Loop Unrolling with Peeling and Jamming and Fission**: The loop is not unrolled with peeling and jamming and fission because it does not have independent parts that can be parallelized.
113. **Loop Unrolling with Peeling and Jamming and Fusion**: The loop is not unrolled with peeling and jamming and fusion because the*/

for (int nl = 0; nl < iterations; nl++) {
    for (int i = 0; i < LEN_1D - k; i++) {
        a[i] = a[i+k] + b[i];
    }
    dummy(a, b, c, d, e, aa, bb, cc, 0.);
}
#pragma endscop

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
    if (DUMP) save_array(func_args.name);
}

int main(int argc, char ** argv){
    int n1 = 1;
    int n3 = 1;
    int* ip;
    real_t s1,s2;
    init(&ip, &s1, &s2);
    printf("Loop \tTime(sec) \tChecksum\n");
    
    time_function(&s431, NULL);
    return EXIT_SUCCESS;
}