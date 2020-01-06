#include <stdlib.h>
#include "GraphBLAS.h"  // in addition to other required C headers                                          

// Compute partial BC metric for a subset of source vertices, s, in graph A
GrB_Info BC_update(GrB_Vector *delta, GrB_Matrix A, GrB_Index *s, GrB_Index nsver)
{
  GrB_Index n; 
  GrB_Matrix_nrows(&n, A);                             // n = # of vertices in graph
  GrB_Vector_new(delta,GrB_FP32,n);                    // Vector<float> delta(n)

  // index and value arrays needed to build numsp
  GrB_Index *i_nsver = (GrB_Index*)malloc(sizeof(GrB_Index)*nsver);
  int32_t   *ones    = (int32_t*)  malloc(sizeof(int32_t)*nsver);
  for(int i=0; i<nsver; ++i) {
    i_nsver[i] = i;
    ones[i] = 1;
  }
  
  // numsp: structure holds the number of shortest paths for each node and starting vertex
  // discovered so far.  Initialized to source vertices:  numsp[s[i],i]=1, i=[0,nsver)
  GrB_Matrix numsp;
  GrB_Matrix_new(&numsp,GrB_INT32,n,nsver);
  GrB_Matrix_build(numsp,s,i_nsver,ones,nsver,GrB_PLUS_INT32);
  free(i_nsver); free(ones);

  // frontier: Holds the current frontier where values are path counts.
  // Initialized to out vertices of each source node in s.
  GrB_Matrix frontier;
  GrB_Matrix_new(&frontier,GrB_INT32,n,nsver);
  GrB_extract(frontier,numsp,GrB_NULL,A,GrB_ALL,n,s,nsver,GrB_DESC_RCT0);

  // sigma: stores frontier information for each level of BFS phase.  The memory
  // for an entry in sigmas is only allocated within the do-while loop if needed.
  // n is an upper bound on diameter.
  GrB_Matrix *sigmas = (GrB_Matrix*)malloc(sizeof(GrB_Matrix)*n);
  
  int32_t   d = 0;                                     // BFS level number
  GrB_Index nvals = 0;                                 // nvals == 0 when BFS phase is complete
  
  // --------------------- The BFS phase (forward sweep) --------------------------- 
  do {
    // sigmas[d](:,s) = d^th level frontier from source vertex s
    GrB_Matrix_new(&(sigmas[d]),GrB_BOOL,n,nsver);
    
    GrB_apply(sigmas[d],GrB_NULL,GrB_NULL,
              GrB_IDENTITY_BOOL,frontier,GrB_NULL);    // sigmas[d](:,:) = (Boolean) frontier 
    GrB_eWiseAdd(numsp,GrB_NULL,GrB_NULL,GrB_PLUS_INT32,
                 numsp,frontier,GrB_NULL);             // numsp += frontier (accum path counts)
    GrB_mxm(frontier,numsp,GrB_NULL,GrB_PLUS_TIMES_SEMIRING_INT32,
            A,frontier,GrB_DESC_RCT0);                 // f<!numsp> = A' +.* f (update frontier)
    GrB_Matrix_nvals(&nvals,frontier);                 // number of nodes in frontier at this level
    d++;
  } while (nvals);

  // nspinv: the inverse of the number of shortest paths for each node and starting vertex.
  GrB_Matrix nspinv;
  GrB_Matrix_new(&nspinv,GrB_FP32,n,nsver);
  GrB_apply(nspinv,GrB_NULL,GrB_NULL,
            GrB_MINV_FP32,numsp,GrB_NULL);             // nspinv = 1./numsp

  // bcu: BC updates for each vertex for each starting vertex in s
  GrB_Matrix bcu;
  GrB_Matrix_new(&bcu,GrB_FP32,n,nsver);
  GrB_assign(bcu,GrB_NULL,GrB_NULL,
             1.0f,GrB_ALL,n,GrB_ALL,nsver,GrB_NULL);  // filled with 1 to avoid sparsity issues

  GrB_Matrix w;                                        // temporary workspace matrix
  GrB_Matrix_new(&w,GrB_FP32,n,nsver);
  
  // -------------------- Tally phase (backward sweep) --------------------
  for (int i=d-1; i>0; i--)  {
    GrB_eWiseMult(w,sigmas[i],GrB_NULL,
                  GrB_TIMES_FP32,bcu,nspinv,GrB_DESC_R);  // w<sigmas[i]>=(1 ./ nsp).*bcu

    // add contributions by successors and mask with that BFS level's frontier
    GrB_mxm(w,sigmas[i-1],GrB_NULL,GrB_PLUS_TIMES_SEMIRING_FP32,
            A,w,GrB_DESC_R);                              // w<sigmas[i-1]> = (A +.* w)
    GrB_eWiseMult(bcu,GrB_NULL,GrB_PLUS_FP32,GrB_TIMES_FP32,
                  w,numsp,GrB_NULL);                      // bcu += w .* numsp
  }
  
  // row reduce bcu and subtract "nsver" from every entry to account 
  // for 1 extra value per bcu row element.
  GrB_reduce(*delta,GrB_NULL,GrB_NULL,GrB_PLUS_FP32,bcu,GrB_NULL);
  GrB_apply(*delta,GrB_NULL,GrB_NULL,GrB_MINUS_FP32,*delta,(float)nsver,GrB_NULL);

  // Release resources
  for(int i=0; i<d; i++) {
    GrB_free(&(sigmas[i]));
  } 
  free(sigmas);
  
  GrB_free(&frontier);     GrB_free(&numsp);
  GrB_free(&nspinv);       GrB_free(&bcu);       GrB_free(&w);
  
  return GrB_SUCCESS;
}
