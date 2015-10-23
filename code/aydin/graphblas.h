#ifndef _GraphBLAS_h_
#define _GraphBLAS_h_

#include <stdint.h>
#include <functional>
#include <vector>
#include <cassert>
#include <algorithm>
#include <iterator>
#include <ostream>
#include <iostream>
#include <cassert>
using namespace std;

/* Copyright © 2015 
 * The Regents of the University of California
 * All Rights Reserved
 * Created by Aydin Buluc, LBNL
 * October 22, 2015
 * abuluc@lbl.gov
 */

// example:
// A = sparse([1 2 3],[1 2 3],[1 1 1],3,3);
// B = sparse([1 2 3],[1 2 3],[1 1 1],3,3);

// C = A * B;
// [i j v] = find(C);



template <typename T>
T CumulativeSum (T * arr, T size)
{
    T prev;
    T tempnz = 0 ;
    for (T i = 0 ; i < size ; ++i)
    {
        prev = arr[i];
        arr[i] = tempnz;
        tempnz += prev ;
    }
    return (tempnz) ;		    // return sum
}

namespace GraphBLAS
{
	template<typename T>
  	class Matrix
  	{
    public:
        Matrix() {};
        Matrix(uint64_t m, uint64_t n):nnz(0), rows(m), cols(n) {};
        Matrix(uint64_t length, uint64_t m, uint64_t n):nnz(length), rows(m), cols(n) {};

        uint64_t rows;
        uint64_t cols;
        uint64_t nnz; // number of nonzeros
        
        uint64_t * colptr;
        uint64_t * rowids;
        T * values;
        
        void PrintVitals()
        {
            cout << "In CSC format: " << endl;
            std::copy(colptr, colptr+ cols+1, ostream_iterator<uint64_t>(cout," ")); cout << endl;
            std::copy(rowids, rowids+ nnz, ostream_iterator<uint64_t>(cout," ")); cout << endl;
            std::copy(values, values+ nnz, ostream_iterator<T>(cout," ")); cout << endl;
        }
        
  	};

    template<typename T>
    Matrix<T> buildMatrix(uint64_t * rowinds, uint64_t * colinds, T * vals, uint64_t length, uint64_t m, uint64_t n)
    {
        // Constructing empty objects (size = 0) are not allowed.
        assert(length != 0 && m != 0 && n != 0);
        
        Matrix<T> M(length, m, n);
        M.colptr = new uint64_t[M.cols+1](); // zero initialize
        M.rowids = new uint64_t[M.nnz];
        M.values = new T[M.nnz];
        
        vector< pair<uint64_t,T> > tosort (M.nnz);
        
        uint64_t * work = new uint64_t[M.cols]();	// workspace
        
        for (uint64_t k = 0 ; k < M.nnz ; ++k)
        {
            uint64_t tmp =  colinds[k];
            work [ tmp ]++ ;		// column counts (i.e, w holds the "col difference array")
        }
        
        if(M.nnz > 0)
        {
            M.colptr[M.cols] = CumulativeSum (work, M.cols) ;		// cumulative sum of w
            copy(work, work+M.cols, M.colptr);
            uint64_t last;
            for (uint64_t k = 0 ; k < M.nnz ; ++k)
            {
                tosort[ work[colinds[k]]++] = make_pair( rowinds[k], vals[k]);
            }
#pragma omp parallel for
            for(uint64_t i=0; i< M.cols; ++i)
            {
                sort(tosort.begin() + M.colptr[i], tosort.begin() + M.colptr[i+1]);
                
                typename vector<pair<uint64_t,T>>::iterator itr;	// iterator is a dependent name
                uint64_t ind;
                for(itr = tosort.begin() + M.colptr[i], ind = M.colptr[i]; itr != tosort.begin() + M.colptr[i+1]; ++itr, ++ind)
                {
                    M.rowids[ind] = itr->first;
                    M.values[ind] = itr->second;
                }
            }
        }
        delete [] work;
        return M;
    }
    
    template<typename T>
    void extractTuples(Matrix<T> & M, uint64_t* & rowinds, uint64_t* & colinds, T* & vals)
    {
        rowinds = new uint64_t[M.nnz];
        colinds = new uint64_t[M.nnz];
        vals = new T[M.nnz];
        
        uint64_t k = 0;
        for(uint64_t i=0; i< M.cols; ++i)
        {
            for(uint64_t j = M.colptr[i]; j < M.colptr[i+1]; ++j)
            {
                rowinds[k] = M.rowids[j];
                colinds[k] = i;
                vals[k++] = M.values[j];
            }
        }
    }
    
    enum Transform      // Input argument preprocessing functions.
    {
        argDesc_null =	0,
        argDesc_neg  =	1,
        argDesc_T 	 =	2,
        argDesc_negT =	3,
        argDesc_notT =	4
    };
    
    // template< class R, class... Args >
    // class function<R(Args...)>
    template <class T1, class T2, class OUT>
    class Semiring
    {
    public:
        OUT SAID; // additive identity
        std::function<OUT(T1, T2)> multiply;
        std::function<OUT(OUT, OUT)> add;
        std::function<OUT(OUT)> alpha;
    };
    
    
    /* C = \alpha A*B
    * SR.alpha does generalized scaling of A*B
    *     SR.alpha = bind2nd(std::multiplies<NT>(), scaler) for constant "scaler", this gives us the BLAS usage
    *     SR.alpha = bind2nd(std::less<NT>(), threshold)?_1:0 for constant threshold, then it works like a drop threshold */
    template <typename OUT, typename T1, typename T2>
    Matrix<OUT> mXm(Semiring<T1, T1, OUT> & SR, Matrix<T1>& A, Matrix<T2>& B)
    {
        assert(A.cols == B.rows);
        Matrix<OUT> C(A.rows, B.cols);
        vector<OUT> spavals(A.rows, SR.SAID);
        vector<bool> spabools(A.rows, false);
        
        vector<uint64_t> * RowIdsofC = new vector<uint64_t>[B.cols];      // row ids for each column of C
        vector<OUT> * ValuesofC = new vector<OUT>[B.cols];      // values for each column of C
        
        for(uint64_t i=0; i < B.cols; ++i)        // for all columns of B
        {
            vector<uint64_t> spainds;
            for(uint64_t j=B.colptr[i]; j < B.colptr[i+1]; ++j)    // For all the nonzeros of the ith column
            {
                uint64_t inner = B.rowids[j];				// get the row id of B (or column id of A)
                for(uint64_t k = A.colptr[inner]; k < A.colptr[inner+1]; ++k)
                {
                    if(!spabools[A.rowids[k]])
                    {
                        spabools[A.rowids[k]] = true;
                        spavals[A.rowids[k]] = SR.multiply(A.values[k], B.values[j]);
                        spainds.push_back(A.rowids[k]);
                    }
                    else
                    {
                        spavals[A.rowids[k]] = SR.add(A.values[k], spavals[A.rowids[k]]);
                    }
                }
            }
            sort(spainds.begin(), spainds.end());
            for(auto itr = spainds.begin(); itr != spainds.end(); ++itr)
            {
                RowIdsofC[i].push_back(*itr);
                ValuesofC[i].push_back(spavals[*itr]);
                spabools[*itr] = false;
            }
        }
        
        C.colptr = new uint64_t[C.cols+1];
        C.colptr[0] = 0;
            
        for(uint64_t i=0; i < C.cols; ++i)        // for all edge lists (do NOT parallelize)
        {
            C.colptr[i+1] = C.colptr[i] + RowIdsofC[i].size();
        }
        C.nnz = C.colptr[C.cols];
        C.rowids = new uint64_t[C.nnz];
        C.values = new OUT[C.nnz];
            
#pragma omp parallel for
        for(uint64_t i=0; i< C.cols; ++i)         // combine step
        {
            transform(ValuesofC[i].begin(), ValuesofC[i].end(), ValuesofC[i].begin(), SR.alpha);
            copy(RowIdsofC[i].begin(), RowIdsofC[i].end(), C.rowids + C.colptr[i]);
            copy(ValuesofC[i].begin(), ValuesofC[i].end(), C.values + C.colptr[i]);
        }
        delete [] RowIdsofC;
        delete [] ValuesofC;
        return C;
    }
}

#endif 
