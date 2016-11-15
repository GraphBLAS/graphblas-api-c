#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include "GraphBLAS.h"

// assign a random number to each element scaled by the inverse of the node's degree.
// This will increase the probability that low degree nodes are selected.
float setRandom(uint32_t degree)
{
  return (0.0001f + Random()/(1.0 + 2.0*degree));   // add 1 to prevent divide by zero
}

/*
 * A variant of Luby's randomized algorithm [Luby 1985]. 
 *
 * Given a numeric n x n adjacency matrix A of an unwieghted and undirected graph (where
 * the value true represents an edge), compute a maximal set of independent vertices and 
 * return it in a boolean n-vector, 'iset' where set[i] == true implies vertex i is a member
 * of the set (the iset vector should be uninitialized on input.)
 */
GrB_info MIS(GrB_Vector *iset, const GrB_Matrix A)
{
  GrB_index n;
  GrB_Matrix_nrows(&n,A);			// n = # of rows of A

  GrB_Vector prob;				// holds random probabilities for each node
  GrB_Vector neighbor_max;			// holds value of max neighbor probability
  GrB_Vector new_members;			// holds set of new members to iset
  GrB_Vector new_neighbors;			// holds set of new neighbors to new iset mbrs.
  GrB_Vector candidates;			// candidate members to iset
  
  GrB_Vector_new(&prob,GrB_FLOAT,n);
  GrB_Vector_new(&neighbor_max,GrB_FLOAT,n);
  GrB_Vector_new(&new_members,GrB_BOOL,n);
  GrB_Vector_new(&new_neighbors,GrB_BOOL,n);
  GrB_Vector_new(&candidates,GrB_BOOL,n);
  GrB_Vector_new(iset,GrB_BOOL,n);		// Initialize independent set vector, bool
  
  GrB_Monoid Max;
  GrB_Monoid_new(&Max,GrB_FLOAT,GrB_MAX_F32, 0.0);
  
  GrB_Semiring maxSelect2nd;                    // Max/Select2nd could this be Max/Times?
  GrB_Semiring_new(&maxSelect2nd, Max, GrB_SECOND_F32);

  GrB_Monoid Lor;
  GrB_Monoid_new(&Lor,GrB_BOOL,GrB_LOR,false);
          
  GrB_Semiring Boolean;				// Boolean semiring
  GrB_Semiring_new(&Boolean,Lor,GrB_LAND);

  // structural complement of mask
  GrB_Descriptor negate_mask;
  GrB_Descriptor_new(&negate_mask);
  GrB_Descriptor_set(&negate_mask,GrB_MASK,GrB_SCMP);

  GrB_UnaryFunction set_random;
  GrB_UnaryFunction_new(&set_random,GrB_UINT32,GrB_FLOAT,setRandom);
  
  // compute the degree of each vertex.
  GrB_Vector degrees;
  GrB_Vector_new(&degrees,GrB_DOUBLE,n);
  GrB_reduce(&degrees,GrB_NULL,GrB_NULL,GrB_PLUS_F64,A,GrB_NULL);

  // candidates[:] = 1.0, (i.e., fill)
  GrB_assign(&candidates,GrB_NULL,GrB_NULL,1.0,GrB_ALL, GrB_NULL); 
  
  // Iterate while there are candidates to check.
  GrB_Index nnz;
  do {
    // compute a random probability scaled by inverse of degree
    GrB_apply(&prob,candidates,GrB_NULL,set_random,degrees,GrB_NULL);
    
    // compute the max probability of all neighbors (could mask = candidates?)
    GrB_mxv(&neighbor_max,GrB_NULL,GrB_NULL,maxSelect2nd,A,prob,GrB_NULL);

    // select vertex if its probability is larger than all its active neighbors
    GrB_eWiseAdd(&new_members,GrB_NULL,GrB_NULL,GrB_GT_DOUBLE,prob,neighbor_max,GrB_NULL);
    
    // add new members to independent set.
    GrB_eWiseAdd(iset,GrB_NULL,GrB_NULL,GrB_LOR,*iset,new_members,GrB_NULL);
    
    // remove new members from set of candidates c = c & !new
    GrB_eWiseMult(&candidates,new_members,GrB_NULL,
                  GrB_LAND,candidates,candidates,negate_mask);
    
    GrB_Vector_nnz(&nnz, candidates);
    if (nnz == 0) { break; } // early exit condition
    
    // Neighbors of new members can also be removed from candidates (mask = candidates?)
    GrB_mxv(&new_neighbors,GrB_NULL,GrB_NULL,Boolean,A,new_members,GrB_NULL);
    GrB_eWiseMult(&candidates,new_neighbors,GrB_NULL,
                  GrB_LAND,candidates,candidates,negate_mask);

    GrB_Vector_nnz(&nnz, candidates);
  } while (nnz);				// if there are no more candidates, we are done.

  GrB_free(neighbor_max);			// free all objects "new'ed"
  GrB_free(new_members);
  GrB_free(new_neighbors);
  GrB_free(prob);
  GrB_free(candidates);
  GrB_free(maxSelect2nd);
  GrB_free(Boolean);
  GrB_free(Max);
  GrB_free(Lor);
  GrB_free(negate_mask);
  GrB_free(set_random);
  GrB_free(degrees);

  return GrB_SUCCESS;
}
