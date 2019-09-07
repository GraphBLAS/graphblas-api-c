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

  // create index ramp for index_of() functionality
  GrB_Index *idx = (GrB_Index*)malloc(N*sizeof(GrB_Index));
  for (GrB_Index i = 0; i < N; ++i) idx[i] = i;
  GrB_Vector index_ramp;
  GrB_Vector_new(&index_ramp, GrB_UINT64, N);
  GrB_Vector_build_UINT64(index_ramp, idx, idx, N, GrB_PLUS_INT64);
  free(idx);

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
    GrB_eWiseMult(wavefront, GrB_NULL, GrB_NULL, GrB_FIRST_UINT64,
                  index_ramp, wavefront, GrB_NULL);

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
  GrB_free(&index_ramp);
   
  return GrB_SUCCESS;
}
