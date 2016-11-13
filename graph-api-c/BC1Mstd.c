#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include "GraphBLAS.h"
#include "GrB_stddef.h"

GrB_info BC(GrB_Vector *delta, GrB_Matrix A, GrB_index s)
/*
 * Given a boolean n x n adjacency matrix A and a source vertex s, 
 * compute the BC-metric vector delta, which should be empty on input.
 */
{
  GrB_index n; 
  GrB_Matrix_nrows(&n, A); 			// n = # of vertices in graph

  GrB_Vector_new(delta,GrB_FP32,n);		// Vector<float> delta(n)

  GrB_Matrix sigma; 				// Matrix<int32_t> sigma(n,n)
  GrB_Matrix_new(&sigma, GrB_INT32, n, n);	// sigma[d,k] = shortest path count to node k at level d

  GrB_Vector q; 
  GrB_Vector_new(&q, GrB_INT32, n);		// Vector<int32_t> q(n) of path counts
  GrB_assign(&q, 1, s);				// q[s] = 1 

  GrB_Vector p; 
  GrB_Vector_new(&p, GrB_INT32, n);		// Vector<int32_t> p(n) shortest path counts so far
  GrB_assign(&p, q);				// p = q

  /*
   * BFS phase
   */
  int32_t d = 0;					// BFS level number
  int32_t sum = 0;					// sum == 0 when BFS phase is complete
  do {
    GrB_assign(&sigma,q,d,GrB_ALL,n);			// sigma[d,:] = q
    GrB_vxm(&q,p,GrB_Int32AddMul,q,A,GrB_CompMask);	// q = # paths to nodes reachable from current level
    GrB_eWiseAdd(&p,GrB_Int32AddMul,p,q);		// accumulate path counts on this level
    GrB_reduce(&sum,GrB_Int32Add,q);			// sum path counts at this level
    d++;
  } while (sum);

  /*
   * BC computation phase
   * (t1,t2,t3,t4) are temporary vectors
   */
  GrB_Vector t1; GrB_Vector_new(&t1,GrB_FP32,n);	
  GrB_Vector t2; GrB_Vector_new(&t2,GrB_FP32,n);
  GrB_Vector t3; GrB_Vector_new(&t3,GrB_FP32,n);
  GrB_Vector t4; GrB_Vector_new(&t4,GrB_FP32,n);
  for(int i=d-1; i>0; i--)
  {
    GrB_assign(&t1,1,GrB_ALL);			// t1 = 1+delta
    GrB_eWiseAdd(&t1,GrB_FP32Add,t1,*delta);
    GrB_assign(&t2,sigma,i,GrB_ALL,n);		// t2 = sigma[i,:]
    GrB_eWiseMult(&t2,GrB_DIV_F32,t1,t2);	// t2 = (1+delta)/sigma[i,:]
    GrB_mxv(&t3,GrB_FP32AddMul,A,t2);		// add contributions made by successors of a node
    GrB_assign(&t4,sigma,i-1,GrB_ALL,n);	// t4 = sigma[i-1,:]
    GrB_eWiseMult(&t4,GrB_FP32Mul,t4,t3);	// t4 = sigma[i-1,:]*t3		
    GrB_eWiseAdd(delta,GrB_FP32Add,*delta,t4);	// accumulate into delta
  }

  GrB_free(&sigma);
  GrB_free(q); GrB_free(p);
  GrB_free(t1); GrB_free(t2); GrB_free(t3); GrB_free(t4);

  return GrB_SUCCESS;
}
