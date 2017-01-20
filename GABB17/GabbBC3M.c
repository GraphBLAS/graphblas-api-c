#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include "GraphBLAS.h"

// s is the list of source vertices
// nsver is the number of source vertices (i.e. the length of s)
GrB_info BC_update(GrB_Vector *delta, GrB_Matrix A, GrB_index *s, GrB_index nsver)
{
  GrB_index n; 
  GrB_Matrix_nrows(&n, A); 			// n = # of vertices in graph

  GrB_Vector_new(delta,GrB_FP32,n);		// Vector<float> delta(n)

  GrB_index *tilln = malloc(sizeof(GrB_index)*nsver);
  GrB_INT32 *ones = malloc(sizeof(GrB_INT32)*nsver);

  for(int i=0; i<nsver; ++i)  {
    tilln[i] = i;
    ones[i] = 1;
  }  // There must be a better way to generate this?
    
  GrB_Matrix frontier; // its nonzero structure holds the current frontier
  GrB_Matrix_new(&frontier, GrB_INT32, n, nsver);	// frontier also stores path counts
  GrB_buildMatrix(&frontier,GrB_NULL,GrB_NULL,s,tilln,ones,nsver,GrB_PLUS_I32,GrB_NULL);  // frontier[s] = 1
    
  GrB_Matrix numsp;     // its nonzero structure holds all vertices that have been discovered so far
  GrB_Matrix_new(&numsp, GrB_INT32, n, nsver);		// numsp also stores shortest path counts so far
  GrB_apply(&numsp,GrB_NULL,GrB_NULL,GrB_IDENTITY_INT32,frontier);	 			// numsp = frontier initially

  GrB_Monoid Int32Add;				// Monoid <int32_t,+,0>
  GrB_Monoid_new(&Int32Add,GrB_INT32,GrB_PLUS_I32,0);
  GrB_Semiring Int32AddMul;			// Semiring <int32_t,int32_t,int32_t,+,*,0,1>
  GrB_Semiring_new(&Int32AddMul,Int32Add,GrB_TIMES_I32);

  GrB_Descriptor desc;				// Descriptor for vxm
  GrB_Descriptor_new(&desc);
  GrB_Descriptor_set(desc,GrB_MASK,GrB_SCMP);	// structural complement of the mask
  GrB_Descriptor_set(desc,GrB_INP0,GrB_TRAN);	// use the transpose of A in mxm below

  GrB_Matrix * sigmas = malloc(sizeof(GrB_Matrix)*n);   // n is an upper bound on diameter
  // The memory for an entry in sigmas is only allocated within the do-while loop if needed
    
  //  BFS phase
  int32_t d = 0;				// BFS level number
  int32_t sum = 0;				// sum == 0 when BFS phase is complete
  do {
    GrB_Matrix_new(&(sigmas[d]), GrB_INT32, n, nsver); // Matrix<int32_t> sigma(n,nsver)
    // sigmas[d](k,s) = shortest path count to node k from starting vertex s at d^th level
      
    GrB_assign(&(sigmas[d]),GrB_NULL,GrB_NULL,frontier,d,GrB_ALL,n,GrB_NULL);// sigma[d,:] = q (aydin: needs update)
      
    GrB_mxm(&q,p,GrB_NULL,Int32AddMul,q,A,desc);		// q = # paths to nodes reachable from current level
    GrB_eWiseAdd(&p,GrB_NULL,GrB_NULL,Int32AddMul,p,q,GrB_NULL);// accumulate path counts on this level
    GrB_reduce(&sum,GrB_NULL,Int32Add,q,GrB_NULL);		// sum path counts at this level
    d++;
  } while (sum);

  // BC computation phase    .... (t1,t2,t3,t4) are temporary vectors
  GrB_Semiring FP32AddMul;	// Semiring <float,float,float,+,*,0.0,1.0>
  GrB_Semiring_new(&FP32AddMul,GrB_FP32,GrB_FP32,GrB_FP32,GrB_PLUS_F32,GrB_TIMES_F32,0.0,1.0);

  GrB_Monoid FP32Add;		// Monoid <float,float,float,+,0.0>
  GrB_Monoid_new(&FP32Add,GrB_FP32,GrB_FP32,GrB_FP32,GrB_PLUS_F32,0.0);

  GrB_Monoid FP32Mul;				// Monoid <float,float,float,*,1.0>
  GrB_Monoid_new(&FP32Mul,GrB_FP32,GrB_FP32,GrB_FP32,GrB_TIMES_F32,1.0);

  GrB_Vector t1; GrB_Vector_new(&t1,GrB_FP32,n);	
  GrB_Vector t2; GrB_Vector_new(&t2,GrB_FP32,n);
  GrB_Vector t3; GrB_Vector_new(&t3,GrB_FP32,n);
  GrB_Vector t4; GrB_Vector_new(&t4,GrB_FP32,n);
  for(int i=d-1; i>0; i--)
  {
    GrB_assign(&t1,GrB_NULL,GrB_NULL,1,GrB_ALL,GrB_NULL);		// t1 = 1+delta
    GrB_eWiseAdd(&t1,GrB_NULL,GrB_NULL,FP32Add,t1,*delta,GrB_NULL);
    GrB_assign(&t2,GrB_NULL,GrB_NULL,sigma,i,GrB_ALL,n,GrB_NULL);	// t2 = sigma[i,:]
    GrB_eWiseMult(&t2,GrB_NULL,GrB_NULL,GrB_DIV_F32,t1,t2,GrB_NULL);	// t2 = (1+delta)/sigma[i,:]
    GrB_mxv(&t3,GrB_NULL,GrB_NULL,FP32AddMul,A,t2,GrB_NULL);		// add contributions made by successors of a node
    GrB_assign(&t4,GrB_NULL,GrB_NULL,sigma,i-1,GrB_ALL,n,GrB_NULL);	// t4 = sigma[i-1,:]
    GrB_eWiseMult(&t4,GrB_NULL,GrB_NULL,FP32Mul,t4,t3,GrB_NULL);	// t4 = sigma[i-1,:]*t3		
    GrB_eWiseAdd(delta,GrB_NULL,GrB_NULL,FP32Add,*delta,t4,GrB_NULL);	// accumulate into delta
  }
  GrB_free(&sigma);
  GrB_free(q); GrB_free(p);
  GrB_free(Int32AddMul); GrB_free(Int32Add); GrB_free(FP32AddMul); 
  GrB_free(FP32Add); GrB_free(FP32Mul);
  GrB_free(desc);
  GrB_free(t1); GrB_free(t2); GrB_free(t3); GrB_free(t4);

  return GrB_SUCCESS;
}