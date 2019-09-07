#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include "GraphBLAS.h"

// Assign a random number to each element scaled by the inverse of the node's degree.
// This will increase the probability that low degree nodes are selected and larger
// sets are selected.
void setRandom(void *out, const void *in)
{ 
  uint32_t degree = *(uint32_t*)in;
  *(float*)out = (0.0001f + random()/(1. + 2.*degree)); // add 1 to prevent divide by zero
}

/*
 * A variant of Luby's randomized algorithm [Luby 1985]. 
 *
 * Given a numeric n x n adjacency matrix A of an unweighted and undirected graph (where
 * the value true represents an edge), compute a maximal set of independent vertices and 
 * return it in a boolean n-vector, 'iset' where set[i] == true implies vertex i is a member
 * of the set (the iset vector should be uninitialized on input.)
 */
GrB_Info MIS(GrB_Vector *iset, const GrB_Matrix A)
{
  GrB_Index n;
  GrB_Matrix_nrows(&n,A);                       // n = # of rows of A

  GrB_Vector prob;                              // holds random probabilities for each node
  GrB_Vector neighbor_max;                      // holds value of max neighbor probability
  GrB_Vector new_members;                       // holds set of new members to iset
  GrB_Vector new_neighbors;                     // holds set of new neighbors to new iset mbrs.
  GrB_Vector candidates;                        // candidate members to iset
  
  GrB_Vector_new(&prob,GrB_FP32,n);
  GrB_Vector_new(&neighbor_max,GrB_FP32,n);
  GrB_Vector_new(&new_members,GrB_BOOL,n);
  GrB_Vector_new(&new_neighbors,GrB_BOOL,n);
  GrB_Vector_new(&candidates,GrB_BOOL,n);
  GrB_Vector_new(iset,GrB_BOOL,n);              // Initialize independent set vector, bool

  GrB_UnaryOp set_random;
  GrB_UnaryOp_new(&set_random,setRandom,GrB_FP32,GrB_UINT32);
  
  // compute the degree of each vertex.
  GrB_Vector degrees;
  GrB_Vector_new(&degrees,GrB_FP64,n);
  GrB_reduce(degrees,GrB_NULL,GrB_NULL,GrB_PLUS_FP64,A,GrB_NULL);

  // Isolated vertices are not candidates: candidates[degrees != 0] = true
  GrB_assign(candidates,degrees,GrB_NULL,true,GrB_ALL,n,GrB_NULL); 

  // add all singletons to iset: iset[degree == 0] = 1
  GrB_assign(*iset,degrees,GrB_NULL,true,GrB_ALL,n,GrB_DESC_RC) ; 

  // Iterate while there are candidates to check.
  GrB_Index nvals;
  GrB_Vector_nvals(&nvals, candidates);
  while (nvals > 0) {
    // compute a random probability scaled by inverse of degree
    GrB_apply(prob,candidates,GrB_NULL,set_random,degrees,GrB_DESC_R);
    
    // compute the max probability of all neighbors
    GrB_mxv(neighbor_max,candidates,GrB_NULL,GrB_MAX_SECOND_SEMIRING_FP32,A,prob,GrB_DESC_R);

    // select vertex if its probability is larger than all its active neighbors,
    // and apply a "masked no-op" to remove stored falses
    GrB_eWiseAdd(new_members,GrB_NULL,GrB_NULL,GrB_GT_FP64,prob,neighbor_max,GrB_NULL);
    GrB_apply(new_members,new_members,GrB_NULL,GrB_IDENTITY_BOOL,new_members,GrB_DESC_R);

    // add new members to independent set.
    GrB_eWiseAdd(*iset,GrB_NULL,GrB_NULL,GrB_LOR,*iset,new_members,GrB_NULL);
    
    // remove new members from set of candidates c = c & !new
    GrB_eWiseMult(candidates,new_members,GrB_NULL,
                  GrB_LAND,candidates,candidates,GrB_DESC_RC);
    
    GrB_Vector_nvals(&nvals, candidates);
    if (nvals == 0) { break; }                  // early exit condition
    
    // Neighbors of new members can also be removed from candidates
    GrB_mxv(new_neighbors,candidates,GrB_NULL,GrB_LOR_LAND_SEMIRING_BOOL,
            A,new_members,GrB_NULL);
    GrB_eWiseMult(candidates,new_neighbors,GrB_NULL,GrB_LAND,
                  candidates,candidates,GrB_DESC_RC);

    GrB_Vector_nvals(&nvals, candidates);
  }

  GrB_free(&neighbor_max);                      // free all objects "new'ed"
  GrB_free(&new_members);
  GrB_free(&new_neighbors);
  GrB_free(&prob);
  GrB_free(&candidates);
  GrB_free(&set_random);
  GrB_free(&degrees);

  return GrB_SUCCESS;
}
