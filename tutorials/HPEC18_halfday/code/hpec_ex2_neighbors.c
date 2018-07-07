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
 * @file hpec_ex2_neighbors.c
 *
 * @brief Code to find neighbors of a vertex; in/out degree
 *
 */

#include <stdio.h>
#include <GraphBLAS.h>
#include "hpec_utils.h"

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
    GrB_Matrix_build(graph, row_indices, col_indices, values, NUM_EDGES,
                     GrB_LOR);

    pretty_print_matrix_UINT64(graph, "GRAPH");

    // Build a wavefront vector with one source node set.
    GrB_Index const SRC_NODE = 6;
    GrB_Vector wavefront;
    GrB_Vector_new(&wavefront, GrB_BOOL, NUM_NODES);
    GrB_Vector_setElement(wavefront, true, SRC_NODE);

    // Build the semiring
    GrB_Monoid Lor;
    GrB_Monoid_new(&Lor, GrB_LOR, (bool)false);
    GrB_Semiring LogicalSR;
    GrB_Semiring_new(&LogicalSR, Lor, GrB_LAND);

    // Build the transpose (INP0) descriptor
    GrB_Descriptor desc_t0;
    GrB_Descriptor_new(&desc_t0);
    GrB_Descriptor_set(desc_t0, GrB_INP0, GrB_TRAN);
    GrB_Descriptor_set(desc_t0, GrB_OUTP, GrB_REPLACE); // not necessary w/o accum

    // find neighbors of SRC_NODE
    // No mask; what if there are self loops?

    pretty_print_vector_UINT64(wavefront, "Source wavefront");
    GrB_mxv(wavefront, GrB_NULL, GrB_NULL,
            LogicalSR, graph, wavefront, desc_t0);
    pretty_print_vector_UINT64(wavefront, "Neighbors");

    // Check results
    {
        bool error_found = false;
        GrB_Index nvals;
        GrB_Vector_nvals(&nvals, wavefront);
        if (nvals != 3)
        {
            fprintf(stderr, "ERROR: wrong number of neighbors (!= 3): %ld\n",
                    nvals);
            error_found = true;
        }

        bool val;

        if (GrB_Vector_extractElement(&val, wavefront, 2UL))
        {
            fprintf(stderr, "ERROR: missing neighbor 2.\n");
            error_found = true;
        }
        if (GrB_Vector_extractElement(&val, wavefront, 3UL))
        {
            fprintf(stderr, "ERROR: missing neighbor 3.\n");
            error_found = true;
        }
        if (GrB_Vector_extractElement(&val, wavefront, 4UL))
        {
            fprintf(stderr, "ERROR: missing neighbor 4.\n");
            error_found = true;
        }

        if (!error_found)
        {
            fprintf(stderr, "neighbor test passed.\n");
        }
    }

    // =======================================================================
    // Extra credit: find out-degree of all vertices
    // =======================================================================

    // Build a degree vector to hold result.
    GrB_Vector out_degree;
    GrB_Vector_new(&out_degree, GrB_UINT64, NUM_NODES);

    GrB_reduce(out_degree, GrB_NULL, GrB_NULL,
               GrB_PLUS_UINT64, graph, GrB_NULL);

    pretty_print_vector_UINT64(out_degree, "out degree");

    // check answer
    {
        bool error_found = false;
        GrB_Index nvals;
        uint64_t  od_answer[] = {2, 2, 1, 2, 1, 1, 3};
        GrB_Vector_nvals(&nvals, out_degree);
        if (nvals != NUM_NODES)
        {
            fprintf(stderr, "ERROR: out degree result missing values. nvals = %ld",
                    nvals);
            error_found = true;
        }

        uint64_t ans;
        for (GrB_Index idx = 0; idx < NUM_NODES; ++idx)
        {
            if (GrB_Vector_extractElement(&ans, out_degree, idx))
            {
                fprintf(stderr, "ERROR: missing value for vertex %ld.\n", idx);
                error_found = true;
            }
            else if (ans != od_answer[idx])
            {
                fprintf(stderr, "ERROR: wrong out degree for vertex %ld, %ld\n",
                        idx, ans);
                error_found = true;
            }
        }

        if (!error_found)
        {
            fprintf(stderr, "out-degree test passed.\n");
        }
    }

    // =======================================================================
    // Extra credit: find in-degree of all vertices
    // =======================================================================

    // Build a degree vector to hold result.
    GrB_Vector in_degree;
    GrB_Vector_new(&in_degree, GrB_UINT64, NUM_NODES);

    GrB_reduce(in_degree, GrB_NULL, GrB_NULL,
               GrB_PLUS_UINT64, graph, desc_t0);

    pretty_print_vector_UINT64(in_degree, "in degree");

    // check answer
    {
        bool error_found = false;
        GrB_Index nvals;
        uint64_t  id_answer[] = {1, 1, 3, 2, 2, 2, 1};
        GrB_Vector_nvals(&nvals, in_degree);
        if (nvals != NUM_NODES)
        {
            fprintf(stderr, "ERROR: out degree result missing values. nvals = %ld",
                    nvals);
            error_found = true;
        }

        uint64_t ans;
        for (GrB_Index idx = 0; idx < NUM_NODES; ++idx)
        {
            if (GrB_Vector_extractElement(&ans, in_degree, idx))
            {
                fprintf(stderr, "ERROR: missing value for vertex %ld.\n", idx);
                error_found = true;
            }
            else if (ans != id_answer[idx])
            {
                fprintf(stderr, "ERROR: wrong in degree for vertex %ld, %ld\n",
                        idx, ans);
                error_found = true;
            }
        }

        if (!error_found)
        {
            fprintf(stderr, "in-degree test passed.\n");
        }
    }

    // Cleanup
    GrB_free(&graph);
    GrB_free(&wavefront);
    GrB_free(&Lor);
    GrB_free(&LogicalSR);
    GrB_free(&desc_t0);
    GrB_free(&out_degree);
    GrB_free(&in_degree);
}
