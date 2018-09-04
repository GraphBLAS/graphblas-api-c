/*
 * Copyright (c) 2018 Carnegie Mellon University.
 * Copyright (c) 2018 Intel Corp.
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
 * @file BuildAdjMatTuple.c
 *
 * @brief Build an adjacency matrix tuple by tuple
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

    // initialize GraphBLAS context
    GrB_init(GrB_BLOCKING);

    GrB_Matrix_new(&graph, GrB_UINT64, NUM_NODES, NUM_NODES);

    GrB_Matrix_setElement(graph, 4, 1, 2);  
    GrB_Matrix_setElement(graph, 4, 2, 1); 
    GrB_Matrix_setElement(graph, 2, 0, 1); 
    GrB_Matrix_setElement(graph, 2, 1, 0); 

    pretty_print_matrix_UINT64(graph, "GRAPH");

    GrB_Index nvals;
    GrB_Matrix_nvals(&nvals, graph);
    assert(nvals == 4);  // verify that we see the number of edges I set

    // Cleanup
    GrB_free(&graph);
    GrB_finalize();
}
