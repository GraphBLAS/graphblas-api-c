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
 * @file  hpec_ex3_level_bfs
 *
 * @brief Implement BFS algorithm that computes vertex level (hops+1 from src)
 *
 */

//#include <stdlib.h>
#include <stdio.h>
//#include <stdint.h>
//#include <stdbool.h>
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
GrB_Info BFS(GrB_Vector *v, GrB_Matrix const A, GrB_Index s)
{
    GrB_Index n;
    GrB_Matrix_nrows(&n,A);                       // n = # of rows of A

    GrB_Vector_new(v,GrB_INT32,n);                // Vector<int32_t> v(n)

    GrB_Vector q;                                 // vertices visited in each level
    GrB_Vector_new(&q,GrB_BOOL,n);                // Vector<bool> q(n)
    GrB_Vector_setElement(q,(bool)true,s);        // q[s] = true, false everywhere else

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
        ++d;                                      // next level (start with 1)
        GrB_assign(*v,q,GrB_NULL,d,GrB_ALL,n,GrB_NULL);   // v[q] = d
        GrB_mxv(q,*v,GrB_NULL,Logical,A,q,desc);  // q[!v] = A' ||.&& q; finds all the
        // unvisited successors from current q
        GrB_Vector_nvals(&nvals, q);              // succ = ||(q)
    } while (nvals);                              // if there is no successor in q, we are done.

    GrB_free(&q);                                 // q vector no longer needed
    GrB_free(&Lor);                               // Logical or monoid no longer needed
    GrB_free(&Logical);                           // Boolean semiring no longer needed
    GrB_free(&desc);                              // descriptor no longer needed

    return GrB_SUCCESS;
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

    // find neighbors of SRC_NODE
    // No mask; what if there are self loops?

    pretty_print_vector_UINT64(levels, "Vertex levels (src is at level 1)");

    // Check results
    {
        bool error_found = false;
        GrB_Index nvals;
        int32_t  lvl_answer[] = {2, 3, 2, 1, 4, 3, 4};
        GrB_Vector_nvals(&nvals, levels);
        if (nvals != NUM_NODES)
        {
            fprintf(stderr, "ERROR: levels result missing values. nvals = %ld",
                    nvals);
            error_found = true;
        }

        // Check against a known result
        int32_t ans;
        for (GrB_Index idx = 0; idx < NUM_NODES; ++idx)
        {
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
            fprintf(stderr, "level test passed.\n");
        }
    }

    // Cleanup
    GrB_free(&graph);
    GrB_free(&levels);
    return 0;
}
