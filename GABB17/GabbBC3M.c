#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include "GraphBLAS.h"

GrB_info BC_update(GrB_Vector *delta, GrB_Matrix A, GrB_index *s, GrB_index nsver)
{
  GrB_index n; 
  GrB_Matrix_nrows(&n, A);                               // n = # of vertices in graph
  GrB_Vector_new(delta,GrB_FP32,n);                      // Vector<float> delta(n)

  GrB_index *tilln = malloc(sizeof(GrB_index)*nsver);
  GrB_INT32 *ones = malloc(sizeof(GrB_INT32)*nsver);
  for(int i=0; i<nsver; ++i) {
    tilln[i] = i;
    ones[i] = 1;
  }
  GrB_Matrix frontier;                                   // nonzero structure holds the current frontier, and
  GrB_Matrix_new(&frontier, GrB_INT32, n, nsver);        // also stores path counts.
  GrB_buildMatrix(&frontier,GrB_NULL,GrB_NULL,s,tilln,ones,nsver,GrB_PLUS_I32,GrB_NULL);  // frontier[s] = 1
    
  GrB_Matrix numsp;                                      // its nonzero structure holds all vertices that have
  GrB_Matrix_new(&numsp, GrB_INT32, n, nsver);           // been discovered and stores shortest path counts so far.
  GrB_apply(&numsp,GrB_NULL,GrB_NULL,GrB_IDENTITY_INT32,frontier);      // (deep copy) numsp = frontier

  GrB_Monoid Int32Add;                                   // Monoid <int32_t,+,0>
  GrB_Monoid_new(&Int32Add,GrB_INT32,GrB_PLUS_I32,0);
  GrB_Semiring Int32AddMul;                              // Semiring <int32_t,int32_t,int32_t,+,*,0,1>
  GrB_Semiring_new(&Int32AddMul,Int32Add,GrB_TIMES_I32);

  GrB_Descriptor desc;                                   // Descriptor for BFS phase mxm
  GrB_Descriptor_new(&desc);
  GrB_Descriptor_set(desc,GrB_MASK,GrB_SCMP);            // structural complement of the mask
  GrB_Descriptor_set(desc,GrB_INP0,GrB_TRAN);            // use the transpose of A in mxm below

  // The memory for an entry in sigmas is only allocated within the do-while loop if needed
  GrB_Matrix * sigmas = malloc(sizeof(GrB_Matrix)*n);    // n is an upper bound on diameter
  int32_t d = 0;                                         // BFS level number
  int32_t nnz = 0;                                       // nnz == 0 when BFS phase is complete
  do {   // ------------------------ BFS phase -----------------------------
    // sigmas[d](k,s) = shortest path count to node k from starting vertex s at d^th level
    GrB_Matrix_new(&(sigmas[d]), GrB_INT32, n, nsver);   // Matrix<int32_t> sigma(n,nsver)
      
    GrB_apply(&(sigmas[d]),GrB_NULL,GrB_NULL,GrB_IDENTITY_INT32,frontier);  // (deep copy) sigma[d,:,:] = frontier
    GrB_mxm(&frontier,numsp,GrB_NULL,Int32AddMul,A,frontier,desc);          // update frontier
    GrB_eWiseAdd(&numsp,GrB_NULL,GrB_NULL,Int32AddMul,numsp,frontier,GrB_NULL);// accumulate path counts
    GrB_Matrix_nnz(&nnz,frontier);                       // number of nodes in frontier at this level
    d++;
  } while (nnz);
    
  GrB_Monoid FP32Add;                                    // Monoid <float,+,0.0>
  GrB_Monoid_new(&FP32Add,GrB_FP32,GrB_PLUS_F32,0.0);
  GrB_Semiring FP32AddMul;                               // Semiring <float,float,float,+,*,0.0,1.0>
  GrB_Semiring_new(&FP32AddMul,FP32Add,GrB_TIMES_F32);
  GrB_Monoid FP32Mul;                                    // Monoid <float,*,1.0>
  GrB_Monoid_new(&FP32Mul,GrB_FP32,GrB_TIMES_F32,1.0);

  GrB_Matrix nspinv;                                     // inverse of the number of shortest paths
  GrB_Matrix_new(&nspinv, GrB_FP32, n, nsver);
  GrB_apply(&nspinv,GrB_NULL,GrB_NULL,GrB_MINV_FP32,numsp);                 // nspinv = 1./nsp
    
  GrB_Matrix bcu;                                        // Betweenness Centrality Updates for each starting vertex in s
  GrB_Matrix_new(&bcu, GrB_FP32, n, nsver);
  GrB_assign(&bcu,GrB_NULL,GrB_NULL,1,GrB_ALL,n, GrB_ALL,nsver,GrB_NULL);   // bcu filled with 1

  GrB_Matrix w; GrB_Matrix_new(&w,GrB_FP32,n, nsver);    // temporary workspace matrix
  for(int i=d-1; i>0; i--)  // ------------------------ BC computation phase ------------------------
  {
      GrB_eWiseMult(&w,GrB_NULL,GrB_NULL,FP32Mul,sigmas[i],nspinv,GrB_NULL);
      GrB_eWiseMult(&w,GrB_NULL,GrB_NULL,FP32Mul,w,bcu,GrB_NULL);           // w = sigmas[i] .* (1 ./ nsp) .* bcu
      // add contributions made by successors of a node
      GrB_mxm(&w,GrB_NULL,GrB_NULL,FP32AddMul,A,w,GrB_NULL);                // w = A +.* w

      GrB_eWiseMult(&w,GrB_NULL,GrB_NULL,FP32Mul,w,sigmas[i-1],GrB_NULL);   // w = w .* sigmas[i-1]
      GrB_eWiseMult(&bcu,GrB_NULL,GrB_PLUS_F32,FP32Mul,w,numsp,GrB_NULL);   // bcu += w .* numsp
  }
  GrB_reduce(&delta,GrB_NULL,GrB_PLUS_F32,GrB_PLUS_F32,bcu,GrB_NULL);
  for(int i=0; i<d; i++)  GrB_free(&(sigmas[i]);         // free sigma matrices
  free(sigmas)
  GrB_free_all(frontier,numsp,nspinv,w,bcu,desc);        // free other matrices and descriptor
  GrB_free_all(Int32AddMul,Int32Add,FP32AddMul, FP32Add, FP32Mul);
  return GrB_SUCCESS;
}
