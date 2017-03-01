#include "GraphBLAS.h"  // in addition to other required C headers                                          |\label{line:include}|

GrB_Info BC_update(GrB_Vector *delta, GrB_Matrix A, GrB_Index *s, GrB_Index nsver) // Compute BC metric     |\label{line:sig}|
{
  GrB_Index n; 
  GrB_Matrix_nrows(&n, A);                               // n = # of vertices in graph
  GrB_Vector_new(delta,GrB_FP32,n);                      // Vector<float> delta(n)                          |\label{line:init_output}|

  GrB_Monoid Int32Add;                                   // Monoid <int32_t,+,0>                            |\label{line:int_add}|
  GrB_Monoid_new(&Int32Add,GrB_INT32,GrB_PLUS_INT32,0);
  GrB_Semiring Int32AddMul;                              // Semiring <int32_t,int32_t,int32_t,+,*,0>        |\label{line:int_arithmetic}|
  GrB_Semiring_new(&Int32AddMul,Int32Add,GrB_TIMES_INT32);

  GrB_Descriptor desc_tsr;                               // Descriptor for BFS phase mxm                    |\label{line:bfs_desc}|
  GrB_Descriptor_new(&desc_tsr);
  GrB_Descriptor_set(desc_tsr,GrB_INP0,GrB_TRAN);        // transpose of the adjacency matrix
  GrB_Descriptor_set(desc_tsr,GrB_MASK,GrB_SCMP);        // structural complement of the mask
  GrB_Descriptor_set(desc_tsr,GrB_OUTP,GrB_REPLACE);     // clear output before result is stored in it.    |\label{line:bfs_desc_end}|

  GrB_Index *i_nsver = malloc(sizeof(GrB_Index)*nsver);  // index and value arrays needed to build numsp    |\label{line:numsp_begin}|
  int32_t   *ones    = malloc(sizeof(int32_t)*nsver);
  for(int i=0; i<nsver; ++i) {
    i_nsver[i] = i;
    ones[i] = 1;
  }
  GrB_Matrix numsp;                                      // Its nonzero structure holds all vertices that have
  GrB_Matrix_new(&numsp, GrB_INT32, n, nsver);           // been discovered and stores number of shortest paths so far.
  GrB_Matrix_build(&numsp,GrB_NULL,GrB_NULL,s,i_nsver,ones,nsver,GrB_PLUS_INT32,GrB_NULL); // numsp[s[i],i]=1, i=[0,nsver)
  free(i_nsver); free(ones);                             //                                                 |\label{line:numsp_end}|

  GrB_Matrix frontier;                                   // Holds the current frontier where values are path counts. |\label{line:frontier_begin}|
  GrB_Matrix_new(&frontier, GrB_INT32, n, nsver);        // Initialized to out vertices of each source node in s.
  GrB_extract(&frontier,numsp,GrB_NULL,A,GrB_ALL,n,s,nsver,desc_tsr);                     //                |\label{line:frontier_end}|

  // The memory for an entry in sigmas is only allocated within the do-while loop if needed
  GrB_Matrix *sigmas = malloc(sizeof(GrB_Matrix)*n);     // n is an upper bound on diameter                 |\label{line:sigma_init}|
  int32_t d = 0;                                         // BFS level number
  int32_t nvals = 0;                                     // nvals == 0 when BFS phase is complete
  do {  // --------------------- The BFS phase (forward sweep) ---------------------------                  |\label{line:dowhile}|
    GrB_Matrix_new(&(sigmas[d]), GrB_BOOL, n, nsver);    // sigmas[d](:,s) = d^th level frontier from source vertex s       |\label{line:sigma_new}|
    GrB_apply(&(sigmas[d]),GrB_NULL,GrB_NULL,GrB_IDENTITY_BOOL,frontier,GrB_NULL); // sigmas[d](:,:) = (Boolean) frontier   |\label{line:sigma_set}|
    GrB_eWiseAdd(&numsp,GrB_NULL,GrB_NULL,Int32Add,numsp,frontier,GrB_NULL);       // numsp += frontier (accum path counts) |\label{line:add_paths}|
    GrB_mxm(&frontier,numsp,GrB_NULL,Int32AddMul,A,frontier,desc_tsr);             // f<!numsp> = A' +.* f (update frontier)|\label{line:mxm1}|
    GrB_Matrix_nvals(&nvals,frontier);                   // number of nodes in frontier at this level
    d++;
  } while (nvals);
    
  GrB_Monoid FP32Add;                                    // Monoid <float,+,0.0>                            |\label{line:fp_arithmetic}|
  GrB_Monoid_new(&FP32Add,GrB_FP32,GrB_PLUS_FP32,0.0f);
  GrB_Monoid FP32Mul;                                    // Monoid <float,*,1.0>
  GrB_Monoid_new(&FP32Mul,GrB_FP32,GrB_TIMES_FP32,1.0f);
  GrB_Semiring FP32AddMul;                               // Semiring <float,float,float,+,*,0.0>
  GrB_Semiring_new(&FP32AddMul,FP32Add,GrB_TIMES_FP32);

  GrB_Matrix nspinv;                                     // inverse of the number of shortest paths         |\label{line:nspinv}|
  GrB_Matrix_new(&nspinv,GrB_FP32,n,nsver);
  GrB_apply(&nspinv,GrB_NULL,GrB_NULL,GrB_MINV_FP32,numsp,GrB_NULL);         // nspinv = 1./numsp
    
  GrB_Matrix bcu;                                        // BC updates for each starting vertex in s        |\label{line:bcu_init}|
  GrB_Matrix_new(&bcu,GrB_FP32,n,nsver);
  GrB_assign(&bcu,GrB_NULL,GrB_NULL,1.0f,GrB_ALL,n, GrB_ALL,nsver,GrB_NULL); // filled with 1 to avoid sparsity issues

  GrB_Descriptor desc_r;                                 // Descriptor for 1st ewisemult in tally           |\label{line:desc}|
  GrB_Descriptor_new(&desc_r);
  GrB_Descriptor_set(desc_r,GrB_OUTP,GrB_REPLACE);       // clear output before result is stored in it.
    
  GrB_Matrix w;                                          // temporary workspace matrix
  GrB_Matrix_new(&w,GrB_FP32,n,nsver);
  for (int i=d-1; i>0; i--)  { // -------------------- Tally phase (backward sweep) --------------------    |\label{line:forloop}|
      GrB_eWiseMult(&w,sigmas[i],GrB_NULL,FP32Mul,bcu,nspinv,desc_r);       // w<sigmas[i]>=(1 ./ nsp).*bcu |\label{line:tallyewm1}|
      
      // add contributions by successors and mask with that BFS level's frontier
      GrB_mxm(&w,sigmas[i-1],GrB_NULL,FP32AddMul,A,w,desc_r);               // w<sigmas[i-1]> = (A +.* w)   |\label{line:mxm2}|
      GrB_eWiseMult(&bcu,GrB_NULL,GrB_PLUS_FP32,FP32Mul,w,numsp,GrB_NULL);  // bcu += w .* numsp            |\label{line:accum_bcu}|
  }
  // subtract "nsver" from every entry in delta (1 extra value per bcu element crept in)
  GrB_assign(delta,GrB_NULL,GrB_NULL,-(float)nsver,GrB_ALL,n,GrB_NULL);     // fill with -nsver   |\label{line:compensate}|
  GrB_reduce(delta,GrB_NULL,GrB_PLUS_FP32,GrB_PLUS_FP32,bcu,GrB_NULL);      // add all updates to -nsver    |\label{line:bcu_reduce}|

  for(int i=0; i<d; i++) { GrB_free(sigmas[i]); } free(sigmas);
  GrB_free_all(frontier,numsp,nspinv,w,bcu,desc_tsr,desc_r);   // macro that expands GrB_free() for each parameter
  GrB_free_all(Int32AddMul,Int32Add,FP32AddMul,FP32Add,FP32Mul);
  return GrB_SUCCESS;
}
