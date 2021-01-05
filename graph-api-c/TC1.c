#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include "GraphBLAS.h"

/*
 * Given an n x n boolean adjacency matrix, A, of an undirected graph, computes
 * the number of triangles in the graph.
 */
uint64_t triangle_count(GrB_Matrix A)
{
  GrB_Index n;
  GrB_Matrix_nrows(&n, A);                        // n = # of vertices
  
  // L: NxN, lower-triangular, bool
  GrB_Matrix L;
  GrB_Matrix_new(&L, GrB_BOOL, n, n);
  GrB_select(L, GrB_NULL, GrB_NULL, GrB_TRIL_UINT64T, A, 0UL, GrB_NULL);

  GrB_Matrix C;
  GrB_Matrix_new(&C, GrB_UINT64, n, n);

  GrB_mxm(C, L, GrB_NULL, GrB_PLUS_TIMES_SEMIRING_UINT64, L, L, GrB_NULL); // C<L> = L +.* L

  uint64_t count;
  GrB_reduce(&count, GrB_NULL, GrB_PLUS_MONOID_UINT64, C, GrB_NULL);     // 1-norm of C

  GrB_free(&C);                      // C matrix no longer needed

  return count;
}
