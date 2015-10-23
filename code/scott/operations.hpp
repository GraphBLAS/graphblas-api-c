/*
 * Copyright (c) 2015 Carnegie Mellon University and The Trustees of Indiana
 * University.
 * All Rights Reserved.
 *
 * THIS SOFTWARE IS PROVIDED "AS IS," WITH NO WARRANTIES WHATSOEVER. CARNEGIE
 * MELLON UNIVERSITY AND THE TRUSTEES OF INDIANA UNIVERSITY EXPRESSLY DISCLAIM
 * TO THE FULLEST EXTENT PERMITTED BY LAW ALL EXPRESS, IMPLIED, AND STATUTORY
 * WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT OF PROPRIETARY RIGHTS.
 *
 * This Program is distributed under a BSD license.  Please see LICENSE file or
 * permission@sei.cmu.edu for more information.  DM-0002659
 */

#ifndef GB_OPERATIONS_HPP
#define GP_OPERATIONS_HPP

#pragma once

#include <cstddef>
#include <vector>

#include <graphblas/math/accum.hpp>
#include <graphblas/math/math.hpp>
#include <graphblas/math/monoid.hpp>
#include <graphblas/semiring/semiring.hpp>

#include <graphblas/detail/config.hpp>

namespace GraphBLAS
{
    /**
     * @brief Perform an element wise binary operation that can be optimized
     *        for "add" semantics (OR short circuit logic).
     *
     * If one element of the Monoid (binary) operation is a structural zero
     * then the result is the other element.  If they are both stored values
     * then the result is the result of the specified Monoid.
     *
     * @tparam <MatrixT>  Models the Matrix concept
     * @tparam <MonoidT>  Models a Monoid (zero and binary function).
     * @tparam <AccumT>   Models a binary function
     *
     * @param[in]  a       The addend (left-hand) operand
     * @param[in]  b       The addend (right-hand) operand
     * @param[out] c       The sum (destination) operand
     * @param[in]  monoid  The element-wise operation to combine the operands
     * @param[in]  accum   The function for how to assign to destination (e.g.
     *                     Accum to "add" to destination, or Assign to
     *                     "replace" the destimation.
     *
     * @throw DimensionException  If a, b, and c do not have the same shape
     */
    template<typename AMatrixT,
             typename BMatrixT,
             typename CMatrixT,
             typename MonoidT =
                 GraphBLAS::math::Plus<typename AMatrixT::ScalarType>,
             typename AccumT =
                 GraphBLAS::math::Assign<typename AMatrixT::ScalarType> >
    inline void eWiseAdd(AMatrixT const &a,
                         BMatrixT const &b,
                         CMatrixT       &c,
                         MonoidT         monoid = MonoidT(),
                         AccumT          accum = AccumT());

    /**
     * @brief Perform an element wise binary operation that can be optimized
     *        for "multiply" semantics (AND short circuit logic).
     *
     * If either element of the Monoid (binary) operation is a structural zero
     * then the result is the structural zero.  If they are both stored values
     * then the result is the result of the specified Monoid.
     *
     * @tparam <MatrixT>  Models the Matrix concept
     * @tparam <MonoidT>  Models a Monoid (zero and binary function).
     * @tparam <AccumT>   Models a binary function
     *
     * @param[in]  a       The multiplicand (left-hand) operand
     * @param[in]  b       The multiplier (right-hand) operand
     * @param[out] c       The product (destination) operand
     * @param[in]  monoid  The element-wise operation to combine the operands
     * @param[in]  accum   The function for how to assign to destination (e.g.
     *                     Accum to "add" to destination, or Assign to
     *                     "replace" the destimation.
     *
     * @throw DimensionException  If a, b, and c do not have the same shape
     */
    template<typename AMatrixT,
             typename BMatrixT,
             typename CMatrixT,
             typename MonoidT =
                 GraphBLAS::math::Times<typename AMatrixT::ScalarType>,
             typename AccumT =
                 GraphBLAS::math::Assign<typename AMatrixT::ScalarType> >
    inline void eWiseMult(AMatrixT const &a,
                          BMatrixT const &b,
                          CMatrixT       &c,
                          MonoidT         monoid = MonoidT(),
                          AccumT          accum = AccumT());

    /**
     * @brief Perform matrix-matrix multiply
     *
     * @tparam <MatrixT>   Models the Matrix concept
     * @tparam <SemiringT> Models a Semiring concept (zero add and mult).
     * @tparam <AccumT>    Models a binary function
     *
     * @param[in]  a  The left-hand matrix
     * @param[in]  b  The right-hand matrix
     * @param[out] c  The result matrix.
     * @param[in]  s  The Semiring to use for the matrix multiplication
     * @param[in]  accum   The function for how to assign to destination (e.g.
     *                     Accum to "add" to destination, or Assign to
     *                     "replace" the destimation.
     *
     * @throw DimensionException  If the matrix dimensions for a, b, and c
     *                            are not consistent for matrix multiply
     */
    template<typename AMatrixT,
             typename BMatrixT,
             typename CMatrixT,
             typename SemiringT =
                 GraphBLAS::ArithmeticSemiring<typename AMatrixT::ScalarType>,
             typename AccumT =
                 GraphBLAS::math::Assign<typename AMatrixT::ScalarType> >
    inline void mXm(AMatrixT const &a,
                    BMatrixT const &b,
                    CMatrixT       &c,
                    SemiringT       s     = SemiringT(),
                    AccumT          accum = AccumT());

    /**
     * @brief Perform row vector-matrix multiply.
     *
     * @tparam <VectorT>   Models the Matrix concept (1 x N)
     * @tparam <MatrixT>   Models the Matrix concept
     * @tparam <SemiringT> Models a Semiring concept (zero add and mult).
     * @tparam <AccumT>    Models a binary function
     *
     * @param[in]  a  The left-hand row vector
     * @param[in]  b  The right-hand matrix
     * @param[out] c  The result row vector.
     * @param[in]  s  The Semiring to use for the matrix multiplication
     * @param[in]  accum   The function for how to assign to destination (e.g.
     *                     Accum to "add" to destination, or Assign to
     *                     "replace" the destimation.
     *
     * @todo Currently the Vector concept is a 1xN Matrix
     *
     * @throw DimensionException  If a or c do not have row vector dimension.
     * @throw DimensionException  If the matrix dimensions for a, b, and c
     *                            are not consistent for matrix multiply
     */
    template<typename AVectorT,
             typename BMatrixT,
             typename CVectorT,
             typename SemiringT =
                 GraphBLAS::ArithmeticSemiring<typename AVectorT::ScalarType>,
             typename AccumT =
                 GraphBLAS::math::Assign<typename AVectorT::ScalarType> >
    inline void vXm(AVectorT const &a,
                    BMatrixT const &b,
                    CVectorT       &c,
                    SemiringT       s     = SemiringT(),
                    AccumT          accum = AccumT());

    /**
     * @brief Perform matrix-column vector multiply.
     *
     * @tparam <MatrixT>   Models the Matrix concept
     * @tparam <VectorT>   Models the Matrix concept (M x 1)
     * @tparam <SemiringT> Models a Semiring concept (zero add and mult).
     * @tparam <AccumT>    Models a binary function
     *
     * @param[in]  a  The left-hand matrix
     * @param[in]  b  The right-hand column vector
     * @param[out] c  The result column vector.
     * @param[in]  s  The Semiring to use for the matrix multiplication
     * @param[in]  accum   The function for how to assign to destination (e.g.
     *                     Accum to "add" to destination, or Assign to
     *                     "replace" the destimation.
     *
     * @todo Currently the Vector concept is a Mx1 Matrix
     *
     * @throw DimensionException  If a or c do not have column vector dimension.
     * @throw DimensionException  If the matrix dimensions for a, b, and c
     *                            are not consistent for matrix multiply
     */
    template<typename AMatrixT,
             typename BVectorT,
             typename CVectorT,
             typename SemiringT =
                 GraphBLAS::ArithmeticSemiring<typename AMatrixT::ScalarType>,
             typename AccumT =
                 GraphBLAS::math::Assign<typename AMatrixT::ScalarType> >
    inline void mXv(AMatrixT const &a,
                    BVectorT const &b,
                    CVectorT       &c,
                    SemiringT       s     = SemiringT(),
                    AccumT          accum = AccumT());

    /**
     * @brief Extract a sub-matrix (sub-graph) from a larger matrix (graph).
     *
     * @tparam <MatrixT>     Models the Matrix concept
     * @tparam <RAIterator>  Models a random access iterator
     * @tparam <AccumT>      Models a binary function
     *
     * @param[in]  a  The matrix to extract from
     * @param[in]  i  Iterator into the ordered set of rows to select from a
     * @param[in]  j  Iterator into the ordered set of columns to select from a
     * @param[out] c  Holds the resulting submatrix
     * @param[in]  accum   The function for how to assign to destination (e.g.
     *                     Accum to "add" to destination, or Assign to
     *                     "replace" the destimation.
     *
     * @note The number of rows in c is taken to equal the number of values
     *       iterated over by i.  The number of columns in c is taken to equal
     *       the number values iterated over by j.
     *
     * @todo Need to throw Dimension exception if attempt to access element
     *       outside the dimensions of a.
     */
    template<typename AMatrixT,
             typename CMatrixT,
             typename RAIteratorI,
             typename RAIteratorJ,
             typename AccumT =
                 GraphBLAS::math::Assign<typename AMatrixT::ScalarType> >
    inline void extract(AMatrixT       const &a,
                        RAIteratorI           i,
                        RAIteratorJ           j,
                        CMatrixT             &c,
                        AccumT                accum = AccumT());

    /**
     * @brief Extract a sub-matrix (sub-graph) from a larger matrix (graph).
     * This wrapper is provided for converting IndexArrayType to iterators
     *
     * @tparam <MatrixT>     Models the Matrix concept
     * @tparam <AccumT>      Models a binary function
     *
     * @param[in]  a  The matrix to extract from
     * @param[in]  i  The ordered set of rows to select from a
     * @param[in]  j  The ordered set of columns to select from a
     * @param[out] c  Holds the resulting submatrix
     * @param[in]  accum   The function for how to assign to destination (e.g.
     *                     Accum to "add" to destination, or Assign to
     *                     "replace" the destimation.
     *
     * @throw DimensionException  If the number of rows in c does not equal the
     *                            number of values in i.
     * @throw DimensionException  If the number of columns in c doesn't equal
     *                            the number of values in j.
     *
     * @todo Need to throw Dimension exception if attempt to access element
     *       outside the dimensions of a.
     */
    template<typename AMatrixT,
             typename CMatrixT,
             typename AccumT =
                 GraphBLAS::math::Assign<typename AMatrixT::ScalarType> >
    inline void extract(AMatrixT       const &a,
                        IndexArrayType const &i,
                        IndexArrayType const &j,
                        CMatrixT             &c,
                        AccumT                accum = AccumT());

    /**
     * @brief Assign a matrix to a set of indices in a larger matrix.
     *
     * @tparam <MatrixT>     Models the Matrix concept
     * @tparam <RAIterator>  Models a random access iterator
     * @tparam <AccumT>      Models a binary function
     *
     * @param[in]  a  The matrix to assign from
     * @param[in]  i  Iterator into the ordered set of rows to assign in c
     * @param[in]  j  Iterator into the ordered set of columns to assign in c
     * @param[out] c  The matrix to assign into a subset of
     * @param[in]  accum   The function for how to assign to destination (e.g.
     *                     Accum to "add" to destination, or Assign to
     *                     "replace" the destimation.
     *
     * @note The number of rows in a is taken to equal the number of values
     *       iterated over by i.  The number of columns in a is taken to equal
     *       the number values iterated over by j.
     *
     * @todo Need to throw Dimension exception if attempt to assign to elements
     *       outside the dimensions of c.
     */
    template<typename AMatrixT,
             typename CMatrixT,
             typename RAIteratorI,
             typename RAIteratorJ,
             typename AccumT =
                 GraphBLAS::math::Assign<typename AMatrixT::ScalarType> >
    inline void assign(AMatrixT const    &a,
                       RAIteratorI        i,
                       RAIteratorJ        j,
                       CMatrixT          &c,
                       AccumT             accum = AccumT());

    /**
     * @brief Assign a matrix to a set of indices in a larger matrix.
     * This wrapper is provided for IndexArrayType to iterators.
     *
     * @tparam <MatrixT>     Models the Matrix concept
     * @tparam <AccumT>      Models a binary function
     *
     * @param[in]  a  The matrix to assign from
     * @param[in]  i  The ordered set of rows to assign in c
     * @param[in]  j  The ordered set of columns to assign in c
     * @param[out] c  The matrix to assign into a subset of
     * @param[in]  accum   The function for how to assign to destination (e.g.
     *                     Accum to "add" to destination, or Assign to
     *                     "replace" the destimation.
     *
     * @throw DimensionException  If the number of rows in a does not equal the
     *                            number of values in i.
     * @throw DimensionException  If the number of columns in a doesn't equal
     *                            the number of values in j.
     *
     * @todo Need to throw Dimension exception if attempt to assign to elements
     *       outside the dimensions of c.
     */
    template<typename AMatrixT,
             typename CMatrixT,
             typename AccumT =
                 GraphBLAS::math::Assign<typename AMatrixT::ScalarType> >
    inline void assign(AMatrixT const       &a,
                       IndexArrayType const &i,
                       IndexArrayType const &j,
                       CMatrixT             &c,
                       AccumT                accum = AccumT());

    /**
     * @brief Apply a unary function to all elements of a matrix.
     *
     * @tparam <MatrixT>         Models the Matrix concept
     * @tparam <UnaryFunctionT>  Models a unary function
     * @tparam <AccumT>          Models a binary function
     *
     * @param[in]  a  The matrix to access from
     * @param[out] c  The matrix to assign into
     * @param[in]  f  The function to apply the result to
     * @param[in]  accum   The function for how to assign to destination (e.g.
     *                     Accum to "add" to destination, or Assign to
     *                     "replace" the destimation.
     *
     * @throw DimensionException  If the sizes of a and c do not match
     */
    template<typename AMatrixT,
             typename CMatrixT,
             typename UnaryFunctionT,
             typename AccumT=
                 GraphBLAS::math::Assign<typename AMatrixT::ScalarType> >
    inline void apply(AMatrixT const &a,
                      CMatrixT       &c,
                      UnaryFunctionT  f,
                      AccumT          accum = AccumT());

    /**
     * @brief Apply a reduction operation to each row of a matrix.
     *
     * @tparam <MatrixT>  Models the Matrix concept
     * @tparam <MonoidT>  Models a binary function
     * @tparam <AccumT>   Models a binary function
     *
     * @param[in]  a  The matrix to perform the row reduction on.
     * @param[out] c  The matrix (column vector) to assign the result to.
     * @param[in]  m  The monoid (binary op) used to reduce the row elements
     * @param[in]  accum   The function for how to assign to destination (e.g.
     *                     Accum to "add" to destination, or Assign to
     *                     "replace" the destimation.
     *
     * @throw DimensionException  If the rows a and c do not match, or the
     *                            the number of columns in c is not one.
     */
    template<typename AMatrixT,
             typename CMatrixT,
             typename MonoidT =
                 GraphBLAS::PlusMonoid<typename AMatrixT::ScalarType>,
             typename AccumT =
                 GraphBLAS::math::Assign<typename AMatrixT::ScalarType> >
    inline void rowReduce(AMatrixT const &a,
                          CMatrixT       &c, // vector?
                          MonoidT         m     = MonoidT(),
                          AccumT          accum = AccumT());

    /**
     * @brief Apply a reduction operation to each column of a matrix.
     *
     * @tparam <MatrixT>  Models the Matrix concept
     * @tparam <MonoidT>  Models a binary function
     * @tparam <AccumT>   Models a binary function
     *
     * @param[in]  a  The matrix to perform the column reduction on.
     * @param[out] c  The matrix (row vector) to assign the result to.
     * @param[in]  m  The monoid (binary op) used to reduce the column elements
     * @param[in]  accum   The function for how to assign to destination (e.g.
     *                     Accum to "add" to destination, or Assign to
     *                     "replace" the destimation.
     *
     * @throw DimensionException  If the columns a and c do not match, or the
     *                            the number of rows in c is not one.
     */
    template<typename AMatrixT,
             typename CMatrixT,
             typename MonoidT =
                 GraphBLAS::PlusMonoid<typename AMatrixT::ScalarType>,
             typename AccumT =
                 GraphBLAS::math::Assign<typename AMatrixT::ScalarType> >
    inline void colReduce(AMatrixT const &a,
                          CMatrixT       &c, // vector?
                          MonoidT         m     = MonoidT(),
                          AccumT          accum = AccumT());

    /**
     * @brief  "Flip" the rows and columns of a matrix
     *
     * @tparam <MatrixT>  Models the Matrix concept
     *
     * @param[in]  a  The matrix to flip
     * @param[out] c  The matrix to assign the result to.
     *
     * @throw DimensionException  If the columns of a != rows of c, or
     *                            if the rows of a != columns of c.
     */
    template<typename AMatrixT,
             typename CMatrixT>
    inline void transpose(AMatrixT const &a,
                          CMatrixT       &c);

    /**
     * @brief Output (row, col, value) tuples from Matrix as three vectors.
     *
     * @note In a departure from the STL standards, we require the
     *       output iterators to be random access (RA) iterators,
     *       although we never read from these iterators, only write.
     *       This allows potentially parallel implementations that
     *       will not be linear in the number of tuples.  The RA
     *       iterators must either be prepared to receive the number
     *       of elements equal to @code get_nnz() or they must include
     *       all necessary memory allocation in their random access logic.
     *
     * @tparam <MatrixT>    model GraphBLAS matrix concept
     * @tparam <RAIterator> model a random access iterator
     *
     * @param[in]  a The matrix with the values to be output
     * @param[out] i The @a row indices output iterator (accepts @code
     *               IndexType values)
     * @param[out] j The @a column indices output iterator (accepts @code
     *               IndexType values)
     * @param[out] v The @a values output iterator (accepts @code
     *               AMatrixT::ScalarType )
     */
    template<typename AMatrixT,
             typename RAIteratorIT,
             typename RAIteratorJT,
             typename RAIteratorVT>
    inline void extractTuples(AMatrixT const &a,
                              RAIteratorIT    i,
                              RAIteratorJT    j,
                              RAIteratorVT    v);
    /**
     * @brief Output (row, col, value) tuples from a Matrix as three vectors.
     * This wrapper is provided for IndexArrayType to iterators.
     *
     * @note This version takes specific container types for outputs.
     *       It is provided to (hopefully) mirror the current
     *       GraphBLAS design discussed by the GraphBLAS committee.
     *       The iterator version subsumes this version.
     *
     * @tparam <MatrixT>    model GraphBLAS matrix concept
     *
     * @param[in]  a The matrix with the values to be output
     * @param[out] i The @a row indices index array
     * @param[out] j The @a column indices index array
     * @param[out] v The @a values vector
     */
    template<typename AMatrixT>
    inline void extractTuples(AMatrixT const                             &a,
                              IndexArrayType                             &i,
                              IndexArrayType                             &j,
                              std::vector<typename AMatrixT::ScalarType> &v);

    /**
     * @brief Populate a Matrix with stored values at specified locations
     *
     * @tparam <MatrixT>      Models the Matrix concept
     * @tparam <RAIteratorT>  Models a random access iterator
     * @tparam <AccumT>       Models a binary function
     *
     * @param[out] m The matrix to assign/accum values to
     * @param[in]  i The iterator over the row indices
     * @param[in]  j The iterator over the col indices
     * @param[in]  v The iterator over the values to store
     * @param[in]  n The number of values to assign.
     * @param[in]  accum   The function for how to assign to destination (e.g.
     *                     Accum to "add" to destination, or Assign to
     *                     "replace" the destimation.
     *
     * @todo Need to add a parameter (functor?) to handle duplicate locations
     *
     * @throw DimensionException  If an element of i or j index outisde the
     *                            size of m.
     */
    template<typename MatrixT,
             typename RAIteratorI,
             typename RAIteratorJ,
             typename RAIteratorV,
             typename AccumT =
                 GraphBLAS::math::Assign<typename MatrixT::ScalarType> >
    inline void buildMatrix(MatrixT     &m,
                            RAIteratorI  i,
                            RAIteratorJ  j,
                            RAIteratorV  v,
                            IndexType    n,
                            AccumT       accum = AccumT());
    /**
     * @brief Populate a Matrix with stored values at specified locations
     * This wrapper is provided for IndexArrayType to iterators.
     *
     * @tparam <MatrixT>      Models the Matrix concept
     * @tparam <AccumT>       Models a binary function
     *
     * @param[out] m The matrix to assign/accum values to
     * @param[in]  i The array of row indices
     * @param[in]  j The array of column indices
     * @param[in]  v The array of values to store
     * @param[in]  accum   The function for how to assign to destination (e.g.
     *                     Accum to "add" to destination, or Assign to
     *                     "replace" the destimation.
     *
     * @todo Need to add a parameter (functor?) to handle duplicate locations
     *
     * @throw DimensionException  If the sizes of i, j, and v are not same.
     * @throw DimensionException  If an element of i or j index outisde the
     *                            size of m.
     */
    template<typename MatrixT,
             typename AccumT =
                 GraphBLAS::math::Assign<typename MatrixT::ScalarType> >
    inline void buildMatrix(MatrixT              &m,
                            IndexArrayType const &i,
                            IndexArrayType const &j,
                            std::vector<typename MatrixT::ScalarType> const &v,
                            AccumT                accum = AccumT());
} // GraphBLAS

#include <graphblas/detail/config.hpp>
#define __GB_SYSTEM_OPERATIONS_HEADER <graphblas/system/__GB_SYSTEM_ROOT/operations.hpp>
#include __GB_SYSTEM_OPERATIONS_HEADER
#undef __GB_SYSTEM_OPERATIONS_HEADER

#endif
