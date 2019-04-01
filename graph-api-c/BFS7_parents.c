#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include "GraphBLAS.h"

/*
 * Given a boolean n x n adjacency matrix A and a source vertex s, performs a BFS traversal
 * of the graph and sets v[i] to the level in which vertex i is visited (v[s] == 1).
 * If i is not reachable from s, then v[i] = 0. (Vector v should be empty on input.)
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

  GrB_Vector_new(parents, GrB_UINT64, n);        // Vector<int32_t> v(n) = 0
  GrB_Vector_setElement(*parents, s, s);         // parents[s] = s

  GrB_Vector wavefront;
  GrB_Vector_new(&wavefront,GrB_UINT64,n);       // Vector<bool> q(n) = false
  GrB_Vector_setElement(wavefront, 1UL, s);      // wavefront[s] = true, false everywhere else

  GrB_Descriptor desc_r;
  GrB_Descriptor_new(&desc_r);
  GrB_Descriptor_set(desc_r, GrB_OUTP, GrB_REPLACE);

  GrB_Descriptor desc_csr;
  GrB_Descriptor_new(&desc_csr);
  GrB_Descriptor_set(desc_csr, GrB_MASK, GrB_SCMP);
  GrB_Descriptor_set(desc_csr, GrB_MASK, GrB_STRUCTURE_ONLY);
  GrB_Descriptor_set(desc_csr, GrB_OUTP, GrB_REPLACE);
  
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

    // Select1st because we are left multiplying wavefront rows
    // Masking out the parent list ensures wavefront values do not
    // overlap values already stored in the parent list
    GrB_vxm(wavefront, *parent_list, GrB_NULL, GxB_MIN_FIRST_UINT64,
            wavefront, graph, desc_csr);

    // We don't need to mask here since we did it in mxm.
    // Merges new parents in current wavefront with existing parents
    // parent_list<!parent_list,merge> += wavefront
    GrB_apply(*parent_list, GrB_NULL, GrB_PLUS_UINT64,
              GrB_IDENTITY_UINT64, wavefront, GrB_NULL);

    GrB_Vector_nvals(&nvals, wavefront);
  }

  GrB_free(); 
  GrB_free();
   
  return GrB_SUCCESS;
}
