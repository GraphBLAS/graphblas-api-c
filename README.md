# GraphBLAS C API

The GraphBLAS C API allows users to implement graph algorithms in the language of linear algebra.
It is part of the larger [GraphBLAS effort](https://graphblas.github.io/), and there are bindings for
the GraphBLAS in many languages (see below).

GraphBLAS provides data structures and algorithms, including matrix and vector data structures,
generalized matrix and vector multiplication, and function objects such as binary operators, monoids,
and semirings that can be used to set the arithmetic used in different operations.

## Example
```C
#include <GraphBLAS.h>

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
```

## GraphBLAS Implementations
* [SuiteSparse GraphBLAS](https://github.com/DrTimothyAldenDavis/GraphBLAS)
* [IBM GraphBLAS](https://github.com/IBM/ibmgraphblas)
