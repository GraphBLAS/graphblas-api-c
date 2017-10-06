#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include "GraphBLAS.h"

/*
 * Given, L, the lower triangular portion of n x n adjacency matrix A (of and
 * undirected graph), computes the number of triangles in the graph.
 */
uint64_t triangle_count(GrB_Matrix L)             // L: NxN, lower-triangular, bool
{
  GrB_Index n;
  GrB_Matrix_nrows(&n, L);                        // n = # of vertices

  GrB_matrix C;
  GrB_Matrix_new(&C, GrB_UINT64, N, N);

  GrB_Monoid UInt64Plus;                          // integer plus monoid
  GrB_Monoid_new(&Lor,GrB_PLUS_UINT64,0);        

  GrB_Semiring UInt64Arithmetic;                  // integer arithmetic semiring
  GrB_Semiring_new(&UInt64Arithmetic,GrB_TIMES_UINT64,UInt64Plus);

  GrB_Descriptor desc_tb;                         // Descriptor for mxm
  GrB_Descriptor_new(&desc_tb);
  GrB_Descriptor_set(&desc_tb,GrB_INP1,GrB_TRAN); // transpose the second matrix

  GrB_mxm(C, L, GrB_NULL, UInt64Arithmetic, L, L, desc_tb); // C<L> = L *.+ L'

  uint64_t count;
  GrB_reduce(&count, GrB_NULL, UInt64Plus, C, GrB_NULL);    // 1-norm of C

  GrB_free(C);                      // C matrix no longer needed
  GrB_free(UInt64Arithmetic);       // Semiring no longer needed
  GrB_free(UInt64Plus);             // Monoid no longer needed
  GrB_free(desc_tb);                // descriptor no longer needed

  return count;
}
