#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include "GraphBLAS.h"

/*
 * Given a boolean n x n adjacency matrix A and a source vertex s, performs a BFS traversal
 * of the graph and sets v[i] to the level in which vertex i is visited (v[s] == 1).
 * If i is not reacheable from s, then v[i] = 0. (Vector v should be empty on input.)
 */
GrB_Info BFS(GrB_Vector *v, GrB_Matrix A, GrB_Index s)
{
  GrB_Index n;
  GrB_Matrix_nrows(&n,A);                       // n = # of rows of A

  GrB_Vector_new(v,GrB_INT32,n);                // Vector<int32_t> v(n)

  GrB_Vector q;                                 // vertices visited in each level
  GrB_Vector_new(&q,GrB_BOOL,n);                // Vector<bool> q(n)
  GrB_Vector_setElement(q,(bool)true,s);        // q[s] = true, false everywhere else

  /*
   * BFS traversal and label the vertices.
   */
  int32_t d = 0;                                // d = level in BFS traversal
  bool succ = false;                            // succ == true when some successor found
  do {
    ++d;                                        // next level (start with 1)
    GrB_assign(*v,q,GrB_NULL,d,GrB_ALL,n,GrB_NULL);   // v[q] = d
    GrB_vxm(q,*v,GrB_NULL,GrB_LOR_LAND_SEMIRING_BOOL,
            q,A,GrB_DESC_RC);                   // q[!v] = q ||.&& A ; finds all the
                                                // unvisited successors from current q
    GrB_reduce(&succ,GrB_NULL,GrB_LOR_MONOID_BOOL,
               q,GrB_NULL);                     // succ = ||(q)
  } while (succ);                               // if there is no successor in q, we are done.

  GrB_free(&q);                                 // q vector no longer needed

  return GrB_SUCCESS;
}
