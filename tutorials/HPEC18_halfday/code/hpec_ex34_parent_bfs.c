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

#include <limits.h>
#include <stdio.h>
#include <GraphBLAS.h>
#include "hpec_utils.h"

//****************************************************************************
void subtract_1(void *out, void const *in)
{
    *(GrB_Index*)out = *(GrB_Index*)in - 1UL;
}

//****************************************************************************
/**
 * @brief Perform a single breadth first traversal on the graph from
 * the specified source vertex.
 *
 * @param[out] parent_list  The list of parents for each traversal (row)
 *                          specified in the roots array. (uninitialized on call)
 * @param[in]  graph        N x N adjacency matrix of the graph on which to
 *                          perform a BFS. (NOT the transpose).  The value
 *                          1 should indicate an edge.
 * @param[in]  source       initial (root) vertex to use in the
 *                          calculation (0-based index).
 * @todo should source be 1-based index?
 */
void bfs_parents(GrB_Vector        *parent_list,
                 GrB_Matrix const   graph,
                 GrB_Index          source)
{
    GrB_Index N;
    GrB_Matrix_nrows(&N, graph);

    GrB_Vector_new(parent_list, GrB_UINT64, N);
    GrB_Vector_setElement(*parent_list, source + 1, source); // parent(src) = src

    GrB_Vector wavefront;
    GrB_Vector_new(&wavefront, GrB_UINT64, N);
    GrB_Vector_setElement(wavefront, (uint64_t)1, source);

    // build a one-based index-of vector
    GrB_Vector index_of;
    GrB_Vector_new(&index_of, GrB_UINT64, N);
    for (GrB_Index ix = 0; ix < N; ++ix)
    {
        GrB_Vector_setElement(index_of, ix+1, ix);
    }

    /// @todo Do we use Suitesparse extension GxB_MIN_SECOND_UINT64 or
    /// do we build the Semiring here?
    GrB_Monoid MinMonoidUINT64;
    GrB_Monoid_new(&MinMonoidUINT64, GrB_MIN_UINT64, (uint64_t)ULONG_MAX);
    GrB_Semiring MinSecondUINT64;
    GrB_Semiring_new(&MinSecondUINT64, MinMonoidUINT64, GrB_SECOND_UINT64);

    GrB_Descriptor desc_r;
    GrB_Descriptor_new(&desc_r);
    GrB_Descriptor_set(desc_r, GrB_OUTP, GrB_REPLACE);

    GrB_Descriptor desc_crt;
    GrB_Descriptor_new(&desc_crt);
    GrB_Descriptor_set(desc_crt, GrB_MASK, GrB_SCMP);
    GrB_Descriptor_set(desc_crt, GrB_INP0, GrB_TRAN);
    GrB_Descriptor_set(desc_crt, GrB_OUTP, GrB_REPLACE);

    GrB_Index nvals;
    GrB_Vector_nvals(&nvals, wavefront);

    while (nvals > 0)
    {
        // convert all stored values to their 1-based index
        GrB_eWiseMult(wavefront, GrB_NULL, GrB_NULL, GrB_SECOND_UINT64,
                      wavefront, index_of, GrB_NULL);

        // Select1st because we are left multiplying wavefront rows
        // Masking out the parent list ensures wavefront values do not
        // overlap values already stored in the parent list
        GrB_mxv(wavefront, *parent_list, GrB_NULL, MinSecondUINT64,
                graph, wavefront, desc_crt);

        // We don't need to mask here since we did it in mxv.
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

    GrB_free(&sub1);
    GrB_free(&desc_crt);
    GrB_free(&desc_r);
    GrB_free(&MinSecondUINT64);
    GrB_free(&MinMonoidUINT64);
    GrB_free(&index_of);
    GrB_free(&wavefront);
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
    GrB_Matrix_build(graph, row_indices, col_indices, &values[0], NUM_EDGES,
                     GrB_LOR);

    pretty_print_matrix_UINT64(graph, "GRAPH");

    // Perform breadth-first traversal from one node and compute parents
    GrB_Index const SRC_NODE = 3;
    GrB_Vector parent_list;
    bfs_parents(&parent_list, graph, SRC_NODE);

    // find neighbors of SRC_NODE
    // No mask; what if there are self loops?

    pretty_print_vector_UINT64(parent_list, "Vertex parents (source has self)");

    // Check results
    {
        bool error_found = false;
        GrB_Index nvals;
        GrB_Index parent_answer[] = {3, 0, 3, 3, 1, 2, 1}; // source slot is self
        GrB_Vector_nvals(&nvals, parent_list);
        if (nvals != NUM_NODES)
        {
            fprintf(stderr, "ERROR: wrong number of values in parent_list. nvals = %ld\n",
                    (long)nvals);
            error_found = true;
        }

        // Check against a known result
        GrB_Index ans;
        for (GrB_Index idx = 0; idx < NUM_NODES; ++idx)
        {
            if (GrB_Vector_extractElement(&ans, parent_list, idx))
            {
                fprintf(stderr, "ERROR: missing parent for vertex %ld.\n", (long)idx);
                error_found = true;
            }
            else if (ans != parent_answer[idx])
            {
                fprintf(stderr, "ERROR: wrong in degree for vertex %ld, %ld\n",
                        (long)idx, (long)ans);
                error_found = true;
            }
        }

        if (!error_found)
        {
            fprintf(stderr, "bfs parent test passed.\n");
        }
    }

    // Cleanup
    GrB_free(&graph);
    GrB_free(&parent_list);
    return 0;
}
