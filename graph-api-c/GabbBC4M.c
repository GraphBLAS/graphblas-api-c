#include "GraphBLAS.h"  // in addition to other required C headers                                          |\label{line:include}|

// Compute partial BC metric for a subset of source vertices, s, in graph A  |\label{line:sig}|
GrB_Info BC_update(GrB_Vector *delta, GrB_Matrix A, GrB_Index *s, GrB_Index nsver)
{
  GrB_Index n; 
  GrB_Matrix_nrows(&n, A);                             // n = # of vertices in graph
  GrB_Vector_new(delta,GrB_FP32,n);                    // Vector<float> delta(n)  |\label{line:init_output}|

  GrB_Monoid Int32Add;                                 // Monoid <int32_t,+,0>   |\label{line:int_add}|
  GrB_Monoid_new(&Int32Add,GrB_INT32,GrB_PLUS_INT32,0);
  GrB_Semiring Int32AddMul;                            // Semiring <int32_t,int32_t,int32_t,+,*,0>  |\label{line:int_arithmetic}|
  GrB_Semiring_new(&Int32AddMul,Int32Add,GrB_TIMES_INT32);

  // Descriptor for BFS phase mxm                    |\label{line:bfs_desc}|
  GrB_Descriptor desc_tsr;
  GrB_Descriptor_new(&desc_tsr);
  GrB_Descriptor_set(desc_tsr,GrB_INP0,GrB_TRAN);      // transpose the adjacency matrix
  GrB_Descriptor_set(desc_tsr,GrB_MASK,GrB_SCMP);      // complement the mask
  GrB_Descriptor_set(desc_tsr,GrB_OUTP,GrB_REPLACE);   // clear output before result is stored

  // index and value arrays needed to build numsp    |\label{line:numsp_begin}|
  GrB_Index *i_nsver = malloc(sizeof(GrB_Index)*nsver);
  int32_t   *ones    = malloc(sizeof(int32_t)*nsver);
  for(int i=0; i<nsver; ++i) {
    i_nsver[i] = i;
    ones[i] = 1;
  }
  
  // numsp: structure holds the number of shortest paths for each node and starting vertex
  // discovered so far.  Initialized to source vertices:  numsp[s[i],i]=1, i=[0,nsver)
  GrB_Matrix numsp;
  GrB_Matrix_new(&numsp, GrB_INT32, n, nsver);
  GrB_Matrix_build(numsp,s,i_nsver,ones,nsver,GrB_PLUS_INT32);
  free(i_nsver); free(ones);    |\label{line:numsp_end}|

  // frontier: Holds the current frontier where values are path counts. |\label{line:frontier_begin}|
  // Initialized to out vertices of each source node in s.
  GrB_Matrix frontier;
  GrB_Matrix_dup(&frontier, numsp);  |\label{line:frontier_end}|

  // sigma: stores frontier information for each level of BFS phase.  The memory
  // for an entry in sigmas is only allocated within the do-while loop if needed
  GrB_Matrix *sigmas = malloc(sizeof(GrB_Matrix)*n);   // n is an upper bound on diameter           |\label{line:sigma_init}|
  
  int32_t d = 0;                                       // BFS level number
  int32_t nvals = 0;                                   // nvals == 0 when BFS phase is complete
  
  // --------------------- The BFS phase (forward sweep) ---------------------------                  |\label{line:dowhile}|
  do {
    // sigmas[d](:,s) = d^th level frontier from source vertex s       |\label{line:sigma_new}|
    GrB_Matrix_new(&(sigmas[d]), GrB_BOOL, n, nsver);
    
    GrB_apply(sigmas[d],GrB_NULL,GrB_NULL,
              GrB_IDENTITY_BOOL,frontier,GrB_NULL);    // sigmas[d](:,:) = (Boolean) frontier   |\label{line:sigma_set}|
    GrB_eWiseAdd(numsp,GrB_NULL,GrB_NULL,
                 Int32Add,numsp,frontier,GrB_NULL);    // numsp += frontier (accum path counts) |\label{line:add_paths}|
    GrB_mxm(frontier,numsp,GrB_NULL,
            Int32AddMul,A,frontier,desc_tsr);          // f<!numsp> = A' +.* f (update frontier)|\label{line:mxm1}|
    GrB_Matrix_nvals(&nvals,frontier);                 // number of nodes in frontier at this level
    d++;
  } while (nvals);
    
  GrB_Monoid FP32Add;                                  // Monoid <float,+,0.0>  |\label{line:fp_arithmetic}|
  GrB_Monoid_new(&FP32Add,GrB_FP32,GrB_PLUS_FP32,0.0f);
  GrB_Monoid FP32Mul;                                  // Monoid <float,*,1.0>
  GrB_Monoid_new(&FP32Mul,GrB_FP32,GrB_TIMES_FP32,1.0f);
  GrB_Semiring FP32AddMul;                             // Semiring <float,float,float,+,*,0.0>
  GrB_Semiring_new(&FP32AddMul,FP32Add,GrB_TIMES_FP32);

  // nspinv: the inverse of the number of shortest paths for each node and starting vertex.  |\label{line:nspinv}|
  GrB_Matrix nspinv;
  GrB_Matrix_new(&nspinv,GrB_FP32,n,nsver);
  GrB_apply(nspinv,GrB_NULL,GrB_NULL,
            GrB_MINV_FP32,numsp,GrB_NULL);             // nspinv = 1./numsp

  // bcu: BC updates for each vertex for each starting vertex in s  |\label{line:bcu_init}|
  GrB_Matrix bcu;
  GrB_Matrix_new(&bcu,GrB_FP32,n,nsver);
  GrB_assign(bcu,GrB_NULL,GrB_NULL,
             1.0f,GrB_ALL,n, GrB_ALL,nsver,GrB_NULL);  // filled with 1 to avoid sparsity issues

  // Descriptor used in the tally phase   |\label{line:desc}|
  GrB_Descriptor desc_r;
  GrB_Descriptor_new(&desc_r);
  GrB_Descriptor_set(desc_r,GrB_OUTP,GrB_REPLACE);     // clear output before result is stored
    
  GrB_Matrix w;                                        // temporary workspace matrix
  GrB_Matrix_new(&w,GrB_FP32,n,nsver);
  
  // -------------------- Tally phase (backward sweep) --------------------    |\label{line:forloop}|
  for (int i=d-1; i>0; i--)  {
    GrB_eWiseMult(w,sigmas[i],GrB_NULL,
                  FP32Mul,bcu,nspinv,desc_r);          // w<sigmas[i]>=(1 ./ nsp).*bcu |\label{line:tallyewm1}|

    // add contributions by successors and mask with that BFS level's frontier
    GrB_mxm(w,sigmas[i-1],GrB_NULL,
            FP32AddMul,A,w,desc_r);                    // w<sigmas[i-1]> = (A +.* w) |\label{line:mxm2}|
    GrB_eWiseMult(bcu,GrB_NULL,GrB_PLUS_FP32,
                  FP32Mul,w,numsp,GrB_NULL);           // bcu += w .* numsp   |\label{line:accum_bcu}|
  }
  
  // subtract "nsver" from every entry in delta (account for 1 extra value per bcu element)
  GrB_assign(*delta,GrB_NULL,GrB_NULL,
             -(float)nsver,GrB_ALL,n,GrB_NULL);        // fill with -nsver   |\label{line:compensate}|
  GrB_reduce(*delta,GrB_NULL,GrB_PLUS_FP32,
             GrB_PLUS_FP32,bcu,GrB_NULL);              // add all updates to -nsver    |\label{line:bcu_reduce}|

  // Release resources
  for(int i=0; i<d; i++) {
    GrB_free(&(sigmas[i]));
  } 
  free(sigmas);
  
  GrB_free(&frontier);     GrB_free(&numsp);
  GrB_free(&nspinv);       GrB_free(&bcu);       GrB_free(&w);
  GrB_free(&desc_tsr);     GrB_free(&desc_r); 
  GrB_free(&Int32AddMul);  GrB_free(&Int32Add);
  GrB_free(&FP32AddMul);   GrB_free(&FP32Add);   GrB_free(&FP32Mul);
  
  return GrB_SUCCESS;
}
