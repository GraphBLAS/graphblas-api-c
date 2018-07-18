/*
 * Copyright (c) 2018 Carnegie Mellon University.
 * All Rights Reserved.
 *
 * THIS SOFTWARE IS PROVIDED "AS IS," WITH NO WARRANTIES WHATSOEVER. CARNEGIE
 * MELLON UNIVERSITY EXPRESSLY DISCLAIMS TO THE FULLEST EXTENT PERMITTED BY LAW
 * ALL EXPRESS, IMPLIED, AND STATUTORY WARRANTIES, INCLUDING, WITHOUT
 * LIMITATION, THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, AND NON-INFRINGEMENT OF PROPRIETARY RIGHTS.
 *
 * This Program is distributed under a BSD license.  Please see LICENSE file or
 * permission@sei.cmu.edu for more information.  DM-0002659
 */

/**
 * @file  hpec_ex4_bfs_parents
 *
 * @brief Find a path from any sink to the src specifed for BFS
 *        algorithm that computes vertex level (hops+1 from src)
 *
 */

#include <stdio.h>
#include <GraphBLAS.h>
#include "hpec_utils.h"

//****************************************************************************
// BFS adapted from Appendix B.1 of C API specification
//****************************************************************************

/*
 * Given a boolean n x n adjacency matrix A and a source vertex s,
 * performs a BFS traversal of the graph and sets v[i] to the level in
 * which vertex i is visited (v[s] == 1).  If i is not reacheable from
 * s, then v[i] = 0. (Vector v should be empty on input.)
 */
GrB_Info BFS(GrB_Vector *v, GrB_Matrix const A, GrB_Index src)
{
  GrB_Index n;
  GrB_Matrix_nrows(&n,A);                       // n = # of rows of A

  GrB_Vector_new(v,GrB_INT32,n);                // Vector<int32_t> v(n)

  GrB_Vector q;                                 // vertices visited in each level
  GrB_Vector_new(&q,GrB_BOOL,n);                // Vector<bool> q(n)
  GrB_Vector_setElement(q,(bool)true,src);      // q[src] = true, false everywhere else

  GrB_Monoid Lor;                               // Logical-or monoid
  GrB_Monoid_new(&Lor,GrB_LOR,(bool)false);

  GrB_Semiring Logical;                         // Logical semiring
  GrB_Semiring_new(&Logical,Lor,GrB_LAND);

  GrB_Descriptor desc;                          // Descriptor for vxm
  GrB_Descriptor_new(&desc);
  GrB_Descriptor_set(desc,GrB_MASK,GrB_SCMP);   // invert the mask
  GrB_Descriptor_set(desc,GrB_INP0,GrB_TRAN);   // transpose the matrix
  GrB_Descriptor_set(desc,GrB_OUTP,GrB_REPLACE);// clear the output before assignment

  /*
   * BFS traversal and label the vertices.
   */
  int32_t d = 0;                                // d = level in BFS traversal
  GrB_Index nvals = 0;                          // nvals > 1 when some successor found
  do {
    ++d;                                        // next level (start with 1)
    GrB_assign(*v,q,GrB_NULL,d,GrB_ALL,n,GrB_NULL);   // v[q] = d
    GrB_mxv(q,*v,GrB_NULL,Logical,A,q,desc);    // q[!v] = A' ||.&& q; finds all the
                                                // unvisited successors from current q
    GrB_Vector_nvals(&nvals, q);                // succ = ||(q)
  } while (nvals);                              // if there is no successor in q, we are done.

  GrB_free(&q);                                 // q vector no longer needed
  GrB_free(&Lor);                               // Logical or monoid no longer needed
  GrB_free(&Logical);                           // Boolean semiring no longer needed
  GrB_free(&desc);                              // descriptor no longer needed

  return GrB_SUCCESS;
}

//****************************************************************************
/*
 * Given a boolean n x n adjacency matrix A, a level vector from BFS and a
 * sink vertex s, finds a path to source vertex.
 * of the graph and sets v[i] to the level in which vertex i is visited (v[s] == 1).
 * If i is not reacheable from s, then v[i] = 0. (Vector v should be empty on input.)
 */
GrB_Info compute_parents(GrB_Vector *parents, GrB_Matrix A, GrB_Vector levels, GrB_Index root)
{
    GrB_Index n;
    GrB_Matrix_nrows(&n, A);                       // n = # vertices

    GrB_Matrix D, Asrc, Adst, Atree;
    GrB_Matrix_new(&D, GrB_UINT64, n, n);
    GrB_Matrix_new(&Asrc, GrB_UINT64, n, n);
    GrB_Matrix_new(&Adst, GrB_UINT64, n, n);
    GrB_Matrix_new(&Atree, GrB_BOOL, n, n);

    GrB_Vector_new(parents, GrB_UINT64, n);

    GrB_Monoid PlusMonoid;                        /// @note use Suitesparse GxB_...
    GrB_Monoid_new(&PlusMonoid, GrB_PLUS_UINT64, (uint64_t)0);

    GrB_Semiring ArithmeticSemiring;              /// @note use Suitesparse GxB_...
    GrB_Semiring_new(&ArithmeticSemiring, PlusMonoid, GrB_TIMES_UINT64);

    GrB_Descriptor desc;                          // Descriptor for annihilation
    GrB_Descriptor_new(&desc);
    GrB_Descriptor_set(desc,GrB_OUTP,GrB_REPLACE);// clear the output before assignment

    // Build a diagonal matrix of levels
    /// @note: This could be done with extractTuples from levels and GrB_Matrix_build
    for (GrB_Index idx = 0; idx < n; ++idx)
    {
        GrB_Index lvl;

        if (GrB_SUCCESS == GrB_Vector_extractElement(&lvl, levels, idx))
        {
            GrB_Matrix_setElement(D, lvl, idx, idx);
        }
    }

    pretty_print_matrix_UINT64(D, "diag(levels)");

    GrB_mxm(Asrc, GrB_NULL, GrB_NULL, ArithmeticSemiring, D, A, GrB_NULL);
    pretty_print_matrix_UINT64(Asrc, "FROM: diag(levels) +.* A");

    GrB_mxm(Adst, GrB_NULL, GrB_NULL, ArithmeticSemiring, A, D, GrB_NULL);
    pretty_print_matrix_UINT64(Adst, "TO:   A +.* diag(levels)");

    // Compute the breadth-first traversal tree (graph)
    GrB_eWiseMult(Atree, GrB_NULL, GrB_NULL, GrB_GT_UINT64, Adst, Asrc, GrB_NULL);
    pretty_print_matrix_BOOL(Atree, "Atree (pre Mask)");

    // COMMON IDIOM
    // Annihilate the stored zeros (IMPORTANT: remember the REPLACE flag)
    GrB_apply(Atree, Atree, GrB_NULL, GrB_IDENTITY_BOOL, Atree, desc);
    pretty_print_matrix_BOOL(Atree, "BFS_tree(A, 3)");

    // Don't need to compute parent of each vertex, MXM with diagonal boolean matrix or N MXV's
    //GrB_Matrix_new(&Parent, GrB_BOOL, n, n);
    //GrB_Matrix_new(&Child, GrB_BOOL, n, n);  // is Identity(n)
    //
    //GrB_mxm(Parent, GrB_NULL, GrB_NULL, ArithmeticSemiring, Atree, Child, GrB_NULL);
    //pretty_print_matrix_BOOL(Parent, "Parents = Atree * I_n");
    //
    //GrB_free(&Parent);
    //GrB_free(&Child);

    // extract the parent information directly from Atree (don't need to do Atree * I_n)
    GrB_Index nvals, nvalstmp, *row_indices, *col_indices;
    bool *values;
    GrB_Matrix_nvals(&nvals, Atree);
    row_indices = (GrB_Index*)malloc(nvals*sizeof(GrB_Index));
    col_indices = (GrB_Index*)malloc(nvals*sizeof(GrB_Index));
    values      =      (bool*)malloc(nvals*sizeof(bool));
    nvalstmp    = nvals;

    GrB_Matrix_extractTuples(row_indices, col_indices, values, &nvalstmp, Atree);
    //assert (nvals == nvalstmp)

    // TRICKY: use the column_indices as values, duplication is benign (choose first or second)
    // root won't have a parent
    GrB_Vector_build(*parents, col_indices, row_indices, nvals, GrB_SECOND_UINT64);
    pretty_print_vector_UINT64(*parents, "PARENTS");

    GrB_free(&desc);
    GrB_free(&ArithmeticSemiring);
    GrB_free(&PlusMonoid);
    GrB_free(&Atree);
    GrB_free(&Adst);
    GrB_free(&Asrc);
    GrB_free(&D);
    return GrB_SUCCESS;
}

//****************************************************************************
void find_path_from(GrB_Index current_id, GrB_Vector const parents)
{
    GrB_Index next_id;
    printf("End node: %ld\n", current_id);
    while (GrB_SUCCESS == GrB_Vector_extractElement(&next_id, parents, current_id))
    {
        printf("Prev. id: %ld\n", next_id);
        current_id = next_id;
    }
}

//****************************************************************************
int main(int argc, char** argv)
{
    GrB_Index const NUM_NODES = 7;
    GrB_Index const NUM_EDGES = 12;
    GrB_Index row_indices[] = {0, 0, 1, 1, 2, 3, 3, 4, 5, 6, 6, 6};
    GrB_Index col_indices[] = {1, 3, 4, 6, 5, 0, 2, 5, 2, 2, 3, 4};
    bool values[] = {true, true, true, true, true, true,
                     true, true, true, true, true, true};
    GrB_Matrix graph;

    GrB_Matrix_new(&graph, GrB_BOOL, NUM_NODES, NUM_NODES);
    GrB_Matrix_build(graph, row_indices, col_indices, (bool*)values, NUM_EDGES,
                     GrB_LOR);

    pretty_print_matrix_UINT64(graph, "GRAPH");

    // Perform breadth-first traversal from one node and compute levels
    GrB_Index const SRC_NODE = 3;
    GrB_Vector levels;
    BFS(&levels, graph, SRC_NODE);

    // Compute parent list from level list
    GrB_Vector parents;
    compute_parents(&parents, graph, levels, SRC_NODE);

    pretty_print_vector_UINT64(levels, "Vertex levels (src is at level 1)");
    pretty_print_vector_UINT64(parents, "Vertex parents (src is empty)");

    find_path_from(4UL, parents);

    // Check results
    {
        bool error_found = false;
        GrB_Index nvals;
        int32_t  lvl_answer[]     = {2, 3, 2, 1, 4, 3, 4};
        GrB_Index parent_answer[] = {3, 0, 3,99, 1, 2, 1};
        GrB_Vector_nvals(&nvals, levels);
        if (nvals != NUM_NODES)
        {
            fprintf(stderr, "ERROR: levels result missing values. nvals = %ld",
                    nvals);
            error_found = true;
        }

        // Check against a known result
        for (GrB_Index idx = 0; idx < NUM_NODES; ++idx)
        {
            int32_t ans;
            if (GrB_Vector_extractElement(&ans, levels, idx))
            {
                fprintf(stderr, "ERROR: missing level value for vertex %ld.\n", idx);
                error_found = true;
            }
            else if (ans != lvl_answer[idx])
            {
                fprintf(stderr, "ERROR: wrong in degree for vertex %ld, %d\n",
                        idx, ans);
                error_found = true;
            }

            // Check parent
            if (ans != 1)
            {
                GrB_Index parent_id;
                if (GrB_Vector_extractElement(&parent_id, parents, idx))
                {
                    fprintf(stderr, "ERROR: missing parent value for vertex %ld.\n", idx);
                    error_found = true;
                }
                else if (parent_id != parent_answer[idx])
                {
                    fprintf(stderr, "ERROR: wrong in degree for vertex %ld, %d\n",
                            idx, ans);
                    error_found = true;
                }
            }
        }

        // Proper directed graph check: using the original (src, dst) edge list.
        // Check level(dst) < level(src) + 2
        for (GrB_Index idx = 0; idx < NUM_EDGES; ++idx)
        {
            GrB_Index lvl_src, lvl_dst;
            if (GrB_Vector_extractElement(&lvl_src, levels, row_indices[idx]))
            {
                fprintf(stderr, "ERROR: missing level value for vertex %ld.\n",
                        row_indices[idx]);
                error_found = true;
            }
            else if (GrB_Vector_extractElement(&lvl_dst, levels, col_indices[idx]))
            {
                fprintf(stderr, "ERROR: missing level value for vertex %ld.\n",
                        col_indices[idx]);
                error_found = true;
            }
            else if ((lvl_dst > lvl_src) && (lvl_dst > (lvl_src + 1)))
            {
                fprintf(stderr, "ERROR: level(dst) too high, (src,dst) = (%ld,%ld), levels = (%ld, %ld)\n",
                        row_indices[idx], col_indices[idx], lvl_src, lvl_dst);
                error_found = true;
            }
        }

        if (!error_found)
        {
            fprintf(stderr, "level and parent tests passed.\n");
        }
    }

    // Cleanup
    GrB_free(&graph);
    GrB_free(&levels);
    GrB_free(&parents);
}
