#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <vector>
 
#include <GraphBLAS/GraphBLAS.hpp>
 
int main(int, char**)
{
	// Syntactic sugar
    typedef double ScalarType;
    GraphBLAS::IndexType const NUM_ROWS = 3;
    GraphBLAS::IndexType const NUM_COLS = 3;

    // Note: size of dimensions required in ctor
    GraphBLAS::Matrix<ScalarType> a(NUM_ROWS, NUM_COLS);
    GraphBLAS::Matrix<ScalarType> b(NUM_ROWS, NUM_COLS);
    GraphBLAS::Matrix<ScalarType> c(NUM_ROWS, NUM_COLS);
 
    // Initialize matrices
    GraphBLAS::IndexArrayType i = {0,  1,  2};  // Note: 0-based indexing
    GraphBLAS::IndexArrayType j = {0,  1,  2};
    std::vector<ScalarType>   v = {1., 1., 1.};
 
    GraphBLAS::buildMatrix(a, i, j, v);
    GraphBLAS::buildMatrix(b, i, j, v);
 
    // Matrix - Matrix multiply
    GraphBLAS::mXm(a, b, c);
 
    // Extract the results: get_nnz() method tells us how big
    GraphBLAS::IndexType nnz = c.get_nnz();
    GraphBLAS::IndexArrayType rows(nnz), cols(nnz);
    std::vector<ScalarType> vals(nnz);
 
    GraphBLAS::extractTuples(c, rows, cols, vals);

	// Check that the answer is correct
    GraphBLAS::IndexArrayType i_res = {0, 1, 2};
    GraphBLAS::IndexArrayType j_res = {0, 1, 2};
    std::vector<ScalarType>   v_res = {1, 1, 1};
 
    bool success = true;
    for (GraphBLAS::IndexType ix = 0; ix < vals.size(); ++ix)
    {
		// Note: no semantics defined about the order of the extracted values,
		// so this in n^2 operation (without sorting)
		bool found = false;
        for (GraphBLAS::IndexType iy = 0; iy < v_res.size(); ++iy)
        {
            if ((i_res[iy] == rows[ix]) && (j_res[iy] == cols[ix]))
			{
				found = true;
				if (v_res[iy] != vals[ix])
				{
					success = false;
				}
				break;
			}
        }
		if (!found)
		{
			success = false;
		}
    }
	
	return (success ? 0 : 1);
}

