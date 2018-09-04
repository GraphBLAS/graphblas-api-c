/*
 * Copyright (c) 2018 Carnegie Mellon University.
 * Copyright (c) 2018 Intel Corporation.
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
 * @file BuildGraph.c
 *
 * @brief Build an adjacency matrix and verify our toolchain works
 *
 */

#include <stdio.h>
#include <assert.h>
#include <GraphBLAS.h>
#include "hpec_utils.h"

//****************************************************************************
int main(int argc, char** argv)
{
    GrB_Index const NUM_NODES = 3;
    GrB_Matrix graph;

    // Initialize a GraphBLAS context
    GrB_init(GrB_BLOCKING);

    GrB_Matrix_new(&graph, GrB_UINT64, NUM_NODES, NUM_NODES);

    GrB_Matrix_setElement(graph, 4, 1, 2);  // set 1,2 element to 4

    pretty_print_matrix_UINT64(graph, "GRAPH");

    GrB_Index nvals;
    GrB_Matrix_nvals(&nvals, graph);
    assert(nvals == 1);

    // Cleanup
    GrB_free(&graph);
    GrB_finalize();
}
