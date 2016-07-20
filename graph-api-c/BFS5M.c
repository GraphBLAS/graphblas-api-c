#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include "GraphBLAS.h"

int32_t d = 0;				// d = level (depth) in BFS traversal
int32_t return_level(bool element)
{
  if (element)                          // The if is unnecessary (for illustation?)
    return d;

  return 0;
}

GrB_info BFS(GrB_Vector *v, const GrB_Matrix A, GrB_index s)
/*
 * Given a boolean n x n adjacency matrix A and a source vertex s, performs a BFS traversal
 * of the graph and sets v[i] to the level in which vertex i is visited (v[s] == 1).
 * If i is not reacheable from s, then v[i] = 0. (Vector v should be empty on input.)
 */
{
  GrB_index n;
  GrB_Matrix_nrows(&n,A);			// n = # of rows of A

  GrB_Vector_new(v,GrB_INT32,n);		// Vector<int32_t> v(n) = 0

  GrB_Vector q;					// vertices visited in each level
  GrB_Vector_new(&q,GrB_BOOL,n);		// Vector<bool> q(n) = false
  GrB_assign(&q,GrB_NULL,true,s);		// q[s] = true, false everywhere else

  GrB_Semiring Boolean;				// Boolean semiring <bool,bool,bool,||,&&,false,true>
  GrB_Semiring_new(&Boolean,GrB_BOOL,GrB_BOOL,GrB_BOOL,GrB_LOR,GrB_LAND,false,true);

  GrB_Descriptor desc;				// Descriptor for vxm
  GrB_Descriptor_new(&desc);
  GrB_Descriptor_set(&desc,GrB_MASK,GrB_SCMP);	// use structural complement of mask

  GrB_UnaryFunction apply_level;
  GrB_UnaryFunction_new(&apply_level,GrB_BOOL,GrB_INT32,return_level);
  
  /*
   * BFS traversal and label the vertices.
   */
  while (Vector_nnz(q)) {			// if there is no successor in q, we are done.
    ++d;					// next level (start with 1)
    // We need to revisit this call to assign. We don't have index_of functionality.
    // [Scott shows alternative using apply in his SIAM AN16 slide deck].
    // GrB_assign(v,GrB_NULL,GrB_NULL,d,index_of(q));	// v[q] = d
    GrB_apply(v,GrB_NULL,GrB_ACCUM,apply_level,q);      // v[q] = d

    GrB_vxm(&q,*v,GrB_NULL,Boolean,q,A,desc);	// q[!v] = q ||.&& A ; finds all the unvisited 
						// successors from current q
  }

  GrB_free(q);					// q vector no longer needed
  GrB_free(Boolean);				// Boolean semiring no longer needed
  GrB_free(desc);				// descriptor no longer needed

  return GrB_SUCCESS;
}
