#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <vector>
 
#include <graphblas/graphblas.hpp>
 
#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE ewise
 
#include <boost/test/included/unit_test.hpp>
 
BOOST_AUTO_TEST_SUITE(mxm)
 
 
BOOST_AUTO_TEST_CASE(mxm)
{
   
    graphblas::LilMatrix<int> a(4,4);
    graphblas::LilMatrix<int> b(4,4);
    graphblas::LilMatrix<int> c(4,4);
 
    // buildmatrix
    graphblas::IndexArrayType i = {    1,   2,   2,    3,    3};
    graphblas::IndexArrayType j = {    0,   1,   2,    1,    3};
    std::vector<int>        v_a = { 4213, 234, 242, 1123, 3342};
    std::vector<int>        v_b = {    2,   2,   2,    2,    2};
 
    graphblas::buildmatrix(a, i, j, v_a);
    graphblas::buildmatrix(b, i, j, v_b);
 
    std::cout << "A = " << std::endl << a << std::endl;
    std::cout << "B = " << std::endl << b << std::endl;
 
    // mxm
    graphblas::mxm(a, b, c);
 
    std::cout << "A * B = " << std::endl << c << std::endl;
 
    //extracttuples:
    graphblas::IndexType nnz = c.get_nnz();
    graphblas::IndexArrayType rows(nnz), cols(nnz);
    std::vector<int> vals(nnz);
 
    graphblas::extracttuples(c, rows, cols, vals);
 
    graphblas::IndexArrayType i_res = {   2,   2,   2,    3,    3,    3};
    graphblas::IndexArrayType j_res = {   0,   1 ,  2,    0,    1,    3};
    std::vector<int>          v_res = { 468, 484, 484, 2246, 6684, 6684};
 
 
    BOOST_CHECK_EQUAL_COLLECTIONS(rows.begin(), rows.end(),
                                  i_res.begin(), i_res.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(cols.begin(), cols.end(),
                                  j_res.begin(), j_res.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(vals.begin(), vals.end(),
                                  v_res.begin(), v_res.end());
}
 
BOOST_AUTO_TEST_SUITE_END()
