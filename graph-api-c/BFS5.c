#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include "GraphBLAS.h"

GrB_info BFS(GrB_Vector *v, GrB_Matrix A, GrB_index s)
/*
 * Given a boolean n x n adjacency matrix A and a source vertex s, performs a BFS traversal
 * of the graph and sets v[i] to the level in which vertex i is visited (v[s] == 1).
 * If i is not reacheable from s, then v[i] = 0. (Vector v should be empty on input.)
 */
{
  GrB_index n;
  GrB_Matrix_nrows(&n,A);			// n = # of rows of A

  GrB_Vector_new(v,GrB_INT32,n);		// Vector<int32_t> v(n)
  GrB_assign(v,0);				// v = 0

  GrB_Vector q;					// vertices visited in each level
  GrB_Vector_new(&q,GrB_BOOL,n);		// Vector<bool> q(n)
  GrB_assign(&q,false);
  GrB_assign(&q,true,s); 			// q[s] = true, false everywhere else

  GrB_Algebra Boolean;				// Boolean algebra <bool,bool,bool,||,&&,false,true>
  GrB_Algebra_new(&Boolean,GrB_BOOL,GrB_BOOL,GrB_BOOL,GrB_LOR,GrB_LAND,false,true);

  GrB_Descriptor desc;				// Descriptor for vxm
  GrB_Descriptor_new(&desc);
  GrB_Descriptor_add(desc,GrB_ARG1,GrB_NOP);	// no operation on the vector
  GrB_Descriptor_add(desc,GrB_ARG2,GrB_NOP);	// no operation on the matrix
  GrB_Descriptor_add(desc,GrB_MASK,GrB_LNOT);	// invert the mask

  /*
   * BFS traversal and label the vertices.
   */
  int32_t d = 1;				// d = level in BFS traversal
  bool succ = false;				// succ == true when some successor found
  do {
    GrB_assign(v,d,q);				// v[q] = d
    GrB_vxm(&q,Boolean,q,A,*v,desc);		// q[!v] = q ||.&& A ; finds all the unvisited 
						// successors from current q
    GrB_reduce(&succ,Boolean,q);		// succ = ||(q)
    d++;					// next level
  } while (succ);				// if there is no successor in q, we are done.

  GrB_free(q);					// q vector no longer needed
  GrB_free(Boolean);				// Boolean semiring no longer needed
  GrB_free(desc);				// descriptor no longer needed

  return GrB_SUCCESS;
}
