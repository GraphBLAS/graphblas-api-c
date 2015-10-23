#include "graphblas.h"
#include <vector>
#include <iostream>
using namespace std;

/* Copyright © 2015
 * The Regents of the University of California
 * All Rights Reserved
 * Created by Aydin Buluc, LBNL
 * October 22, 2015
 * abuluc@lbl.gov
 */

int main()
{
	vector<uint64_t> Arows = {0,1,2};
	vector<uint64_t> Acols = {0,1,2};
	vector<double> Avals = {1.0, 1.0, 1.0};


    GraphBLAS::Matrix<double> A = GraphBLAS::buildMatrix(Arows.data(), Acols.data(), Avals.data(), 3, 3, 3);
	GraphBLAS::Matrix<double> B = GraphBLAS::buildMatrix(Arows.data(), Acols.data(), Avals.data(), 3, 3, 3);

    A.PrintVitals();
    
    GraphBLAS::Semiring<double, double, double> RFDD;   // real field double x double
    RFDD.multiply = std::multiplies<double>();
    RFDD.add = std::plus<double>();
    RFDD.SAID = 0.0;
    RFDD.alpha = std::bind(std::multiplies<double>(), std::placeholders::_1, 1.0); // _1 is the first argument passed to the expression
    
    GraphBLAS::Matrix<double> C = GraphBLAS::mXm(RFDD, A, B);

    // have the user allocate memory; that would make ExtractTuples symmetic
    uint64_t * Crows;
    uint64_t * Ccols;
    double * Cvals;
    GraphBLAS::extractTuples(C, Crows, Ccols, Cvals);
    
    cout << "Output in triples:" << endl;
    std::copy(Crows, Crows+ C.nnz, ostream_iterator<uint64_t>(cout," ")); cout << endl;
    std::copy(Ccols, Ccols+ C.nnz, ostream_iterator<uint64_t>(cout," ")); cout << endl;
    std::copy(Cvals, Cvals+ C.nnz, ostream_iterator<double>(cout," ")); cout << endl;

	return 0;
}
