#include "graphblas.h"
#include <vector>
#include <iostream>
using namespace std;

/* Aydin Buluc, LBNL
 * October 2, 2015
 * abuluc@lbl.gov
 */

int main()
{
	vector<uint64_t> Arows = {1,1,1};
	vector<uint64_t> Acols = {2,3,1};
	vector<double> Avals = {3.0, 3.0, 1.0};


    GraphBLAS::Matrix<double> A = GraphBLAS::BuildMatrix(Arows.data(), Acols.data(), Avals.data(), 3, 3, 3);
	GraphBLAS::Matrix<double> B = GraphBLAS::BuildMatrix(Arows.data(), Acols.data(), Avals.data(), 3, 3, 3);

    
    GraphBLAS::Matrix<double> C;    // if dimensions are set during construction, MxM checks them before execution. otherwise MxM sets them.
    auto alphaop = bind2nd(std::plus<double>(), 0.0);   // learn how to wrap it
    GraphBLAS::MxM(C, A, B, multiplies<double>(), plus<double>(), alphaop, 0.0);    // accumulation is missing

    // have the user allocate memory; that would make ExtractTuples symmetic
    uint64_t * Crows;
    uint64_t * Ccols;
    double * Cvals;
    GraphBLAS::ExtractTuples(C, Crows, Ccols, Cvals);
    
    std::copy(Crows, Crows+ C.nnz, ostream_iterator<uint64_t>(cout," ")); cout << endl;
    std::copy(Ccols, Ccols+ C.nnz, ostream_iterator<uint64_t>(cout," ")); cout << endl;
    std::copy(Cvals, Cvals+ C.nnz, ostream_iterator<double>(cout," ")); cout << endl;

	return 0;
}
