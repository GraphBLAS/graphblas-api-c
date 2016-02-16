#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include "GraphBLAS.h"

GrB_info BC(GrB_Vector *delta, GrB_Matrix A, GrB_index s)
/*
 * Given a boolean n x n adjacency matrix A and a source vertex s, 
 * compute the BC-metric vector delta, which should be empty on input.
 */
{
  GrB_index n; 
  GrB_Matrix_nrows(&n, A); 			// n = # of vertices in graph

  GrB_Vector_new(delta,GrB_FP32,n);		// Vector<float> delta(n)
  GrB_assign(delta, 0); 			// delta = 0

  GrB_Matrix sigma; 				// Matrix<int32_t> sigma(n,n)
  GrB_Matrix_new(&sigma, GrB_INT32, n, n);	// sigma[d,k] = shortest path count to node k
  GrB_assign(&sigma, 0);			// at level d

  GrB_Vector q; 
  GrB_Vector_new(&q, GrB_INT32, n);		// Vector<int32_t> q(n) of path counts
  GrB_assign(&q, 0);				// q = 0
  GrB_assign(&q, 1, s);				// q[s] = 1 

  GrB_Vector p; 
  GrB_Vector_new(&p, GrB_INT32, n);		// Vector<int32_t> p(n) shortest path counts so far
  GrB_assign(&p, q);				// p = q

  GrB_Space Int32;				// Space <int32_t,int32_t,int32_t,+,*,0,1>
  GrB_Space_new(&Int32,GrB_INT32,GrB_INT32,GrB_INT32,GrB_PLUS,GrB_TIMES,0,1);

  GrB_Descriptor desc;                          // Descriptor for vxm
  GrB_Descriptor_new(&desc);
  GrB_Descriptor_add(desc,GrB_ARG1,GrB_NOP);    // no operation on the vector
  GrB_Descriptor_add(desc,GrB_ARG2,GrB_NOP);    // no operation on the matrix
  GrB_Descriptor_add(desc,GrB_MASK,GrB_LNOT);   // invert the mask

  /*
   * BFS phase
   */
  int32_t d = 0;				// BFS level number
  int32_t sum = 0;				// sum == 0 when BFS phase is complete
  do {
    GrB_assign(&sigma,q,d,GrB_ALL);		// sigma[d,:] = q
    GrB_vxm(&q,Int32,q,A,p,desc);		// q = # paths to nodes reachable from current level
    GrB_ewiseadd(&p,Int32,p,q);			// accumulate path counts on this level
    GrB_reduce(&sum,q,GrB_PLUS);		// sum path counts at this level
    d++;
  } while (sum);

  /*
   * BC computation phase
   * (t1,t2,t3,t4) are temporary vectors
   */
  GrB_Space FP32AddMul;				// Space <float,float,float,+,*,0,1>
  GrB_Space_new(&FP32AddMul,GrB_FP32,GrB_FP32,GrB_FP32,GrB_PLUS,GrB_TIMES,0.0,1.0);

  GrB_Space FP32AddDiv;				// Space <float,float,float,+,/,0,1>
  GrB_Space_new(&FP32AddDiv,GrB_FP32,GrB_FP32,GrB_FP32,GrB_PLUS,GrB_DIV,0.0,1.0);

  GrB_Vector t1; GrB_Vector_new(&t1,GrB_FP32,n);	
  GrB_Vector t2; GrB_Vector_new(&t2,GrB_FP32,n);
  GrB_Vector t3; GrB_Vector_new(&t3,GrB_FP32,n);
  GrB_Vector t4; GrB_Vector_new(&t4,GrB_FP32,n);
  for(int i=d-1; i>0; i--)
  {
    GrB_assign(&t1, 1);				// t1 = 1+delta
    GrB_ewiseadd(&t1,FP32AddMul,t1,*delta);
    GrB_assign(&t2,sigma,i,GrB_ALL);		// t2 = sigma[i,:]
    GrB_ewisemul(&t2,FP32AddDiv,t1,t2);		// t2 = (1+delta)/sigma[i,:]
    GrB_mxv(&t3,FP32AddMul,A,t2);		// add contributions made by successors of a node
    GrB_assign(&t4,sigma,i-1,GrB_ALL);		// t4 = sigma[i-1,:]
    GrB_ewisemul(&t4,FP32AddMul,t4,t3);		// t4 = sigma[i-1,:]*t3		
    GrB_ewiseadd(delta,FP32AddMul,*delta,t4);	// accumulate into delta
  }

  GrB_free(q); GrB_free(p);
  GrB_free(Int32); GrB_free(FP32AddMul); GrB_free(FP32AddDiv);
  GrB_free(desc);

  return GrB_SUCCESS;
}
