#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include "GraphBLAS.h"

/*
 * Given a binary n x n adjacency matrix A and a source vertex s, performs a BFS
 * traversal of the graph and sets parents[i] to the index of vertex i's parent.
 * The parent of the root vertex, s, will be set to itself (parents[s] == s). If
 * vertex i is not reachable from s, parents[i] will not contain a stored value. 
 */
GrB_Info BFS(GrB_Vector *parents, const GrB_Matrix A, GrB_Index s) 
{
  GrB_Index N;
  GrB_Matrix_nrows(&N, A);                       // N = # vertices

  GrB_Vector_new(parents, GrB_UINT64, N);
  GrB_Vector_setElement(*parents, s, s);         // parents[s] = s

  GrB_Vector wavefront;
  GrB_Vector_new(&wavefront, GrB_UINT64, N);
  GrB_Vector_setElement(wavefront, 1UL, s);      // wavefront[s] = 1

  /*
   * BFS traversal and label the vertices.
   */
  GrB_Index nvals;
  GrB_Vector_nvals(&nvals, wavefront);

  while (nvals > 0)
  {
    // convert all stored values in wavefront to their 0-based index
    GrB_apply(wavefront, GrB_NULL, GrB_NULL, GrB_ROWINDEX_UINT64,
              wavefront, 0UL, GrB_NULL);

    // "FIRST" because left-multiplying wavefront rows. Masking out the parent
    // list ensures wavefront values do not overwrite parents already stored.
    GrB_vxm(wavefront, *parents, GrB_NULL, GrB_MIN_FIRST_SEMIRING_UINT64,
            wavefront, A, GrB_DESC_RSC);

    // Don't need to mask here since we did it in mxm. Merges new parents in
    // current wavefront with existing parents: parents += wavefront
    GrB_apply(*parents, GrB_NULL, GrB_PLUS_UINT64,
              GrB_IDENTITY_UINT64, wavefront, GrB_NULL);

    GrB_Vector_nvals(&nvals, wavefront);
  }

  GrB_free(&wavefront);
   
  return GrB_SUCCESS;
}
