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
 * @file bfs.c
 *
 * @brief Standard breadth first search implementation.
 *
 */

#include <stdio.h>
#include <GraphBLAS.h>
#include "algorithm_utils.h"

//****************************************************************************
void subtract_1(void *out, void const *in)
{
    *(GrB_Index*)out = *(GrB_Index*)in - 1UL;
}

//****************************************************************************
/**
 * @brief Perform a single breadth first search (BFS) traversal on the
 *        given graph.
 *
 * @param[in]  graph        N x N adjacency matrix of the graph on which to
 *                          perform a BFS. (NOT the transpose).  The value
 *                          1 should indicate an edge.
 * @param[in]  roots        N-vector, initial root(s) to use in the
 *                          calculation. It should have a
 *                          single value set to '1' corresponding to the
 *                          root.
 * @param[out] parent_list  The list of parents for each traversal (row)
 *                          specified in the roots array. (uninitialized on call)
 */
void bfs(GrB_Matrix const   graph,
         GrB_Vector const   roots,
         GrB_Vector        *parent_list)
{
    GrB_Vector wavefront;
    GrB_Vector_dup(&wavefront, roots);

    // Set the roots parents to themselves using one-based indices because
    // the mask is sensitive to stored zeros.
    GrB_Vector_dup(parent_list, roots);
    index_of(*parent_list, 1);

    GrB_Descriptor desc_r;
    GrB_Descriptor_new(&desc_r);
    GrB_Descriptor_set(desc_r, GrB_OUTP, GrB_REPLACE);

    GrB_Descriptor desc_cr;
    GrB_Descriptor_new(&desc_cr);
    GrB_Descriptor_set(desc_cr, GrB_MASK, GrB_SCMP);
    GrB_Descriptor_set(desc_cr, GrB_OUTP, GrB_REPLACE);

    GrB_Index nvals;
    GrB_Vector_nvals(&nvals, wavefront);

    while (nvals > 0)
    {
        // convert all stored values to their 1-based index
        index_of(wavefront, 1);

        // Select1st because we are left multiplying wavefront rows
        // Masking out the parent list ensures wavefront values do not
        // overlap values already stored in the parent list
        GrB_vxm(wavefront, *parent_list, GrB_NULL, GxB_MIN_FIRST_UINT64,
                wavefront, graph, desc_cr);

        // We don't need to mask here since we did it in mxm.
        // Merges new parents in current wavefront with existing parents
        // parent_list<!parent_list,merge> += wavefront
        GrB_apply(*parent_list, GrB_NULL, GrB_PLUS_UINT64,
                  GrB_IDENTITY_UINT64, wavefront, GrB_NULL);

        GrB_Vector_nvals(&nvals, wavefront);
    }

    // Restore zero-based indices by subtracting 1 from all stored values
    GrB_UnaryOp sub1;
    GrB_UnaryOp_new(&sub1, subtract_1, GrB_UINT64, GrB_UINT64);

    GrB_apply(*parent_list, GrB_NULL, GrB_NULL, sub1, *parent_list, desc_r);
}

//************************************************************************
/**
 * @brief Perform a set of breadth first search (BFS) traversals on the
 *        given graph.
 *
 * @param[in]  graph        N x N adjacency matrix of the graph on which to
 *                          perform a BFS. (NOT the transpose).  The value
 *                          1 should indicate an edge,
 * @param[in]  roots        R x N, initial root(s) to use in the
 *                          calculation. Each row (1) corresponds to a
 *                          different BFS traversal and (2) should have a
 *                          single value set to '1' corresponding to the
 *                          root.
 * @param[out] parent_list  The list of parents for each traversal (row)
 *                          specified in the roots array. (uninitialized on call)
 */
void bfs_batch(GrB_Matrix const   graph,
               GrB_Matrix const   roots,
               GrB_Matrix        *parent_list)
{
    GrB_Matrix wavefronts;
    GrB_Matrix_dup(&wavefronts, roots);

    // Set the roots parents to themselves using one-based indices because
    // the mask is sensitive to stored zeros.
    GrB_Matrix_dup(parent_list, roots);
    col_index_of(*parent_list, 1);

    GrB_Descriptor desc_r;
    GrB_Descriptor_new(&desc_r);
    GrB_Descriptor_set(desc_r, GrB_OUTP, GrB_REPLACE);

    GrB_Descriptor desc_cr;
    GrB_Descriptor_new(&desc_cr);
    GrB_Descriptor_set(desc_cr, GrB_MASK, GrB_SCMP);
    GrB_Descriptor_set(desc_cr, GrB_OUTP, GrB_REPLACE);

    GrB_Index nvals;
    GrB_Matrix_nvals(&nvals, wavefronts);

    while (nvals > 0)
    {
        // convert all stored values to their 1-based column index
        col_index_of(wavefronts, 1);

        // Select1st because we are left multiplying wavefront rows
        // Masking out the parent list ensures wavefronts values do not
        // overlap values already stored in the parent list
        GrB_mxm(wavefronts, *parent_list, GrB_NULL, GxB_MIN_FIRST_UINT64,
                wavefronts, graph, desc_cr);

        // We don't need to mask here since we did it in mxm.
        // Merges new parents in current wavefront with existing parents
        // parent_list<!parent_list,merge> += wavefronts
        GrB_apply(*parent_list, GrB_NULL, GrB_PLUS_UINT64,
                  GrB_IDENTITY_UINT64, wavefronts, GrB_NULL);

        GrB_Matrix_nvals(&nvals, wavefronts);
    }

    // Restore zero-based indices by subtracting 1 from all stored values
    GrB_UnaryOp sub1;
    GrB_UnaryOp_new(&sub1, subtract_1, GrB_UINT64, GrB_UINT64);
    GrB_apply(*parent_list, GrB_NULL, GrB_NULL, sub1, *parent_list, desc_r);
}

//****************************************************************************
GrB_Index depth = 0;
void return_depth(void *out, void const *in)
{
    *(GrB_Index*)out = depth;
}

//************************************************************************
/**
 * @brief Perform a single breadth first searches (BFS) on the given graph.
 *
 * @param[in]  graph      NxN adjacency matrix of the graph on which to
 *                        perform a BFS (not the transpose).  A value of
 *                        1 indicates an edge (structural zero = 0).
 * @param[in]  roots      N-vector initial roots to use in the calculation
 *                        of R simultaneous traversals.  A value of 1 in a
 *                        given position indicates a root for the
 *                        traversal..
 * @param[out] levels     The level (distance in unweighted graphs) from
 *                        the corresponding root of that BFS.  Roots are
 *                        assigned a value of 1. (a value of 0 implies not
 *                        reachable.
 */
GrB_Info bfs_level_masked(GrB_Matrix const  graph,
                          GrB_Vector const  roots,     //row vector
                          GrB_Vector       *levels)
{
    GrB_Vector wavefront;
    GrB_Vector_dup(&wavefront, roots);

    /// Assert graph is square/have a compatible shape with wavefront
    GrB_Index grows, gcols;
    GrB_Matrix_nrows(&grows, graph);
    GrB_Matrix_ncols(&gcols, graph);

    GrB_Index wsize;
    GrB_Vector_size(&wsize, wavefront);

    if ((grows != gcols) || (wsize != grows))
    {
        return GrB_DIMENSION_MISMATCH;
    }

    GrB_Vector_new(levels, GrB_UINT64, wsize);

    depth = 0;  // note: not thread safe

    GrB_UnaryOp apply_depth;
    GrB_UnaryOp_new(&apply_depth, return_depth, GrB_UINT64, GrB_BOOL);

    GrB_Descriptor desc_cr;
    GrB_Descriptor_new(&desc_cr);
    GrB_Descriptor_set(desc_cr, GrB_OUTP, GrB_REPLACE);
    GrB_Descriptor_set(desc_cr, GrB_MASK, GrB_SCMP);

    GrB_Index nvals;
    GrB_Vector_nvals(&nvals, wavefront);

    while (nvals > 0)
    {
        // Increment the level
        ++depth;

        // Apply the level to all newly visited nodes
        GrB_apply(*levels, GrB_NULL, GrB_PLUS_UINT64,
                  apply_depth, wavefront, GrB_NULL);

        // Advance the wavefront and mask out nodes already assigned levels
        GrB_vxm(wavefront, *levels, GrB_NULL, GxB_LOR_LAND_BOOL,
                wavefront, graph, desc_cr);

        GrB_Vector_nvals(&nvals, wavefront);
    }

    GrB_free(&desc_cr);
    GrB_free(&apply_depth);

    return GrB_SUCCESS;
}

//****************************************************************************
/**
 * @brief Perform a breadth first search (BFS) on the given graph.
 *
 * @param[in]  graph        The graph to perform a BFS on.  NOT built from
 *                          the transpose of the adjacency matrix.
 *                          (1 indicates edge, structural zero = 0).
 * @param[in]  roots        The initial wavefront to use in the calculation.
 *                          (1 indicates root, structural zero = 0).
 * @param[out] levels       The level (distance in unweighted graphs) from
 *                          the corresponding root of that BFS. (uninitialized on call)
 */
void bfs_level_batch_masked(GrB_Matrix const   graph,
                            GrB_Matrix const   roots,       // row vectors
                            GrB_Matrix        *levels)
{
    GrB_Matrix wavefronts;
    GrB_Matrix_dup(&wavefronts, roots);

    /// @todo Assert graph is square
    /// @todo Assert graph has a compatible shape with wavefront?
    GrB_Index rows, cols;
    GrB_Matrix_nrows(&rows, wavefronts);
    GrB_Matrix_ncols(&cols, wavefronts);

    GrB_Matrix_new(levels, GrB_UINT64, rows, cols);

    depth = 0;  // note: not thread safe

    GrB_UnaryOp apply_depth;
    GrB_UnaryOp_new(&apply_depth, return_depth, GrB_UINT64, GrB_BOOL);

    GrB_Descriptor desc_cr;
    GrB_Descriptor_new(&desc_cr);
    GrB_Descriptor_set(desc_cr, GrB_OUTP, GrB_REPLACE);
    GrB_Descriptor_set(desc_cr, GrB_MASK, GrB_SCMP);

    GrB_Index nvals;
    GrB_Matrix_nvals(&nvals, wavefronts);

    while (nvals > 0)
    {
        // Increment the level
        ++depth;

        GrB_apply(*levels, GrB_NULL, GrB_PLUS_UINT64,
                  apply_depth, wavefronts, GrB_NULL);

        GrB_mxm(wavefronts, *levels, GrB_NULL, GxB_LOR_LAND_BOOL,
                wavefronts, graph, desc_cr);

        GrB_Matrix_nvals(&nvals, wavefronts);
    }

    GrB_free(&desc_cr);
    GrB_free(&apply_depth);
}

//****************************************************************************
int main(int argc, char** argv)
{
    // FA Tech Note graph - note that vertex 7 is isolated.
    //
    //    {-, -, -, 1, -, -, -, -, -},
    //    {-, -, -, 1, -, -, 1, -, -},
    //    {-, -, -, -, 1, 1, 1, -, 1},
    //    {1, 1, -, -, 1, -, 1, -, -},
    //    {-, -, 1, 1, -, -, -, -, 1},
    //    {-, -, 1, -, -, -, -, -, -},
    //    {-, 1, 1, 1, -, -, -, -, -},
    //    {-, -, -, -, -, -, -, -, -},
    //    {-, -, 1, -, 1, -, -, -, -};

    GrB_Index const NUM_NODES = 9;
    GrB_Index const NUM_EDGES = 20;
    GrB_Index i[] = {0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3,
                     4, 4, 4, 5, 6, 6, 6, 8, 8};
    GrB_Index j[] = {3, 3, 6, 4, 5, 6, 8, 0, 1, 4, 6,
                     2, 3, 8, 2, 1, 2, 3, 2, 4};
    GrB_Index v[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                     1, 1, 1, 1, 1, 1, 1, 1, 1};

    GrB_Matrix graph;
    GrB_Matrix_new(&graph, GrB_UINT64, NUM_NODES, NUM_NODES);
    GrB_Matrix_build(graph, i, j, v, NUM_EDGES, GrB_PLUS_UINT64);
    pretty_print_matrix_UINT64(graph, "GRAPH");

    // ***********************************************************************
    // run BFS/single root algorithm
    {
        GrB_Vector parent_list;
        GrB_Vector root;
        GrB_Vector_new(&root, GrB_UINT64, NUM_NODES);
        GrB_Vector_setElement(root, 1, 0);

        bfs(graph, root, &parent_list);
        pretty_print_vector_UINT64(parent_list, "parents(GRAPH, root=0)");

        GrB_free(&parent_list);
        GrB_free(&root);
    }

    // ***********************************************************************
    // run BFS/multi-root algorithm
    {
        GrB_Matrix parent_lists;
        GrB_Matrix roots;
        GrB_Matrix_new(&roots, GrB_UINT64, NUM_NODES, NUM_NODES);
        for (GrB_Index ix = 0; ix < NUM_NODES; ++ix)
            GrB_Matrix_setElement(roots, 1, ix, ix);

        bfs_batch(graph, roots, &parent_lists);
        pretty_print_matrix_UINT64(parent_lists, "parents(GRAPH, all roots)");

        GrB_free(&parent_lists);
        GrB_free(&roots);
    }

    // ***********************************************************************
    // run BFS-level/single-root algorithm (using boolean adj matrix)
    {
        GrB_Index const START_INDEX = 5;
        bool vb[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                     1, 1, 1, 1, 1, 1, 1, 1, 1};
        GrB_Matrix bool_graph;
        GrB_Matrix_new(&bool_graph, GrB_BOOL, NUM_NODES, NUM_NODES);
        GrB_Matrix_build(bool_graph, i, j, vb, NUM_EDGES, GrB_LOR);

        GrB_Vector levels;
        GrB_Vector root;
        GrB_Vector_new(&root, GrB_BOOL, NUM_NODES);
        GrB_Vector_setElement(root, true, START_INDEX);

        bfs_level_masked(bool_graph, root, &levels);
        pretty_print_vector_UINT64(levels, "levels(GRAPH, root=5)");

        GrB_free(&levels);
        GrB_free(&root);
        GrB_free(&bool_graph);
    }

    // ***********************************************************************
    // run BFS-level/multi-root algorithm (using boolean adj matrix)
    {
        bool vb[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                     1, 1, 1, 1, 1, 1, 1, 1, 1};
        GrB_Matrix bool_graph;
        GrB_Matrix_new(&bool_graph, GrB_BOOL, NUM_NODES, NUM_NODES);
        GrB_Matrix_build(bool_graph, i, j, vb, NUM_EDGES, GrB_LOR);

        GrB_Matrix levels;
        GrB_Matrix roots;
        GrB_Matrix_new(&roots, GrB_BOOL, NUM_NODES, NUM_NODES);
        for (GrB_Index ix = 0; ix < NUM_NODES; ++ix)
            GrB_Matrix_setElement(roots, true, ix, ix);

        bfs_level_batch_masked(bool_graph, roots, &levels);
        pretty_print_matrix_UINT64(levels, "levels(GRAPH, all roots)");

        GrB_free(&levels);
        GrB_free(&roots);
        GrB_free(&bool_graph);
    }
}
