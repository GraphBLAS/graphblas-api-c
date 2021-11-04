#include <stdlib.h>
#include "GraphBLAS.h"  // in addition to other required C headers                                          |\label{line:include}|

// Compute partial BC metric for a subset of source vertices, s, in graph A  |\label{line:sig}|
GrB_Info BC_update(GrB_Vector *delta, GrB_Matrix A, GrB_Index *s, GrB_Index nsver)
{
  GrB_Index n;
  GrB_Matrix_nrows(&n, A);                             // n = # of vertices in graph
  GrB_Vector_new(delta,GrB_FP32,n);                    // Vector<float> delta(n)  |\label{line:init_output}|

  // Monoid:   GrB_PLUS_MONOID_INT32         <int32_t,+,0>   |\label{line:int_add}|
  // Semiring: GrB_PLUS_TIMES_SEMIRING_INT32 <int32_t,int32_t,int32_t,+,*,0>  |\label{line:int_arithmetic}|

  // Descriptor for BFS phase mxm: GrB_DESC_RST0             |\label{line:bfs_desc}|
  //   transpose the adjacency matrix
  //   complement the mask
  //   clear output before result is stored

  // index and value arrays needed to build numsp    |\label{line:numsp_begin}|
  GrB_Index *i_nsver = (GrB_Index*)malloc(sizeof(GrB_Index)*nsver);
  int32_t   *ones    = (int32_t*)  malloc(sizeof(int32_t)*nsver);
  for(int i=0; i<nsver; ++i) {
    i_nsver[i] = i;
    ones[i] = 1;
  }

  // numsp: structure holds the number of shortest paths for each node and starting vertex
  // discovered so far.  Initialized to source vertices:  numsp[s[i],i]=1, i=[0,nsver)
  GrB_Matrix numsp;
  GrB_Matrix_new(&numsp, GrB_INT32, n, nsver);
  GrB_Matrix_build(numsp,s,i_nsver,ones,nsver,GrB_PLUS_INT32);
  free(i_nsver); free(ones);    			// |\label{line:numsp_end}|

  // frontier: Holds the current frontier where values are path counts. |\label{line:frontier_begin}|
  // Initialized to out vertices of each source node in s.
  GrB_Matrix frontier;
  GrB_Matrix_new(&frontier, GrB_INT32, n, nsver);
  GrB_extract(frontier,numsp,GrB_NULL,A,GrB_ALL,n,s,nsver,GrB_DESC_RST0); // |\label{line:frontier_end}|

  // sigma: stores frontier information for each level of BFS phase.  The memory
  // for an entry in sigmas is only allocated within the do-while loop if needed
  GrB_Matrix *sigmas = (GrB_Matrix*)malloc(sizeof(GrB_Matrix)*n);   // n is an upper bound on diameter           |\label{line:sigma_init}|

  int32_t   d = 0;                                       // BFS level number
  GrB_Index nvals = 0;                                   // nvals == 0 when BFS phase is complete

  // --------------------- The BFS phase (forward sweep) ---------------------------                  |\label{line:dowhile}|
  do {
    // sigmas[d](:,s) = d^th level frontier from source vertex s       |\label{line:sigma_new}|
    GrB_Matrix_new(&(sigmas[d]), GrB_BOOL, n, nsver);

    GrB_apply(sigmas[d],GrB_NULL,GrB_NULL,
              GrB_IDENTITY_BOOL,frontier,GrB_NULL);    // sigmas[d](:,:) = (Boolean) frontier   |\label{line:sigma_set}|
    GrB_eWiseAdd(numsp,GrB_NULL,GrB_NULL,GrB_PLUS_INT32,
                 numsp,frontier,GrB_NULL);             // numsp += frontier (accum path counts) |\label{line:add_paths}|
    GrB_mxm(frontier,numsp,GrB_NULL,GrB_PLUS_TIMES_SEMIRING_INT32,
            A,frontier,GrB_DESC_RST0);                 // f<!numsp> = A' +.* f (update frontier)|\label{line:mxm1}|
    GrB_Matrix_nvals(&nvals,frontier);                 // number of nodes in frontier at this level
    d++;
  } while (nvals);

  // GrB_PLUS_MONOID_FP32:  Monoid <float,+,0.0>  |\label{line:fp_arithmetic}|
  // GrB_TIMES_MONOID_FP32: Monoid <float,*,1.0>
  // GrB_PLUS_TIMES_SEMIRING_FP32: Semiring <float,float,float,+,*,0.0>

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

  // GrB_DESC_R: Descriptor used in the tally phase   |\label{line:desc}|

  GrB_Matrix w;                                        // temporary workspace matrix
  GrB_Matrix_new(&w,GrB_FP32,n,nsver);

  // -------------------- Tally phase (backward sweep) --------------------    |\label{line:forloop}|
  for (int i=d-1; i>0; i--)  {
    GrB_eWiseMult(w,sigmas[i],GrB_NULL,GrB_TIMES_FP32,
                  bcu,nspinv,GrB_DESC_R);              // w<sigmas[i]>=(1 ./ nsp).*bcu |\label{line:tallyewm1}|

    // add contributions by successors and mask with that BFS level's frontier
    GrB_mxm(w,sigmas[i-1],GrB_NULL,GrB_PLUS_TIMES_SEMIRING_FP32,
            A,w,GrB_DESC_R);                           // w<sigmas[i-1]> = (A +.* w) |\label{line:mxm2}|
    GrB_eWiseMult(bcu,GrB_NULL,GrB_PLUS_FP32,
                  GrB_TIMES_FP32,w,numsp,GrB_NULL);    // bcu += w .* numsp   |\label{line:accum_bcu}|
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

  return GrB_SUCCESS;
}
