#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include "GraphBLAS.h"

int32_t level = 0;			// level = depth in BFS traversal, roots=1, unvisited=0
int32_t return_level(bool element)
{
  if (element)                          // The if is unnecessary (for illustation?)
    return level;

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

  GrB_Monoid Lor;
  GrB_Monoid_new(&Lor,GrB_BOOL,GrB_LOR,false);
          
  GrB_Semiring Boolean;				// Boolean semiring
  GrB_Semiring_new(&Boolean,Lor,GrB_LAND);

  GrB_Descriptor scmp_mask;			// Descriptor for structural complement of mask
  GrB_Descriptor_new(&scmp_mask);
  GrB_Descriptor_set(&scmp_masl,GrB_MASK,GrB_SCMP);

  GrB_UnaryOp apply_level;
  GrB_UnaryOp_new(&apply_level,GrB_BOOL,GrB_INT32,return_level);
  
  /*
   * BFS traversal and label the vertices.
   */
  level = 0;
  while (Vector_nnz(q)) {			// if there is no successor in q, we are done.
    ++level;					// next level (start with 1)
    GrB_apply(v,GrB_NULL,GrB_PLUS_I32,apply_level,q,GrB_NULL);  // v[q] = level
    GrB_vxm(&q,*v,GrB_NULL,Boolean,q,A,scmp_mask,GrB_NULL);	// q[!v] = q ||.&& A ; finds all the 
                                                        	// unvisited successors from current q
  }

  GrB_free(q);					// q vector no longer needed
  GrB_free(Lor)
  GrB_free(Boolean);				// Boolean semiring no longer needed
  GrB_free(scmp_mask);				// descriptor no longer needed

  return GrB_SUCCESS;
}