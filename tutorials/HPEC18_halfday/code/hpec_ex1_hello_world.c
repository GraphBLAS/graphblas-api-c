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
 * @file  hpec_ex1_hello_world.c
 *
 * @brief Build a GraphBLAS adjacency matrix for the logo graph and
 *        make sure it is correct.
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

    // Discuss: 1-based vs. 0-based indexing
    GrB_Index row_indices[] = {0, 0, 1, 1, 2, 3, 3, 4, 5, 6, 6, 6};
    GrB_Index col_indices[] = {1, 3, 4, 6, 5, 0, 2, 5, 2, 2, 3, 4};
    bool values[] = {true, true, true, true, true, true,
                     true, true, true, true, true, true};
    GrB_Matrix graph;

    // Discuss options for Domain
    GrB_Matrix_new(&graph, GrB_BOOL, NUM_NODES, NUM_NODES);

    // Discuss options for dup BinaryOp
    GrB_Matrix_build(graph, row_indices, col_indices, values, NUM_EDGES,
                     GrB_LOR);

    // Provided in hpec_utils.h
    pretty_print_matrix_UINT64(graph, "GRAPH");

    // Check results
    bool error_found = false;

    GrB_Index nvals;
    GrB_Matrix_nvals(&nvals, graph);
    printf("Number of edges: %ld\n", nvals);

    if (nvals != NUM_EDGES)
    {
        fprintf(stderr, "ERROR: wrong number of edges (!= %ld): %ld\n",
                NUM_EDGES, nvals);
        error_found = true;
    }
    for (size_t idx = 0; idx < NUM_EDGES; ++idx)
    {
        bool val;

        if (GrB_Matrix_extractElement(&val, graph,
                                      row_indices[idx], col_indices[idx]))
        {
            fprintf(stderr, "ERROR: missing element at: (%ld, %ld)\n",
                    row_indices[idx], col_indices[idx]);
            error_found = true;
        }
        else
        {
            if (val != values[idx])
            {
                fprintf(stderr, "ERROR: wrong value stored at: (%ld, %ld)\n",
                        row_indices[idx], col_indices[idx]);
                error_found = true;
            }
        }
    }

    if (!error_found)
    {
        fprintf(stderr, "All tests passed.\n");
    }

    // Cleanup
    GrB_free(&graph);
}
