//=================================================================================================
/*!
//  \file blaze/math/expressions/TSMatSMatKronExpr.h
//  \brief Header file for the transpose sparse matrix/sparse matrix Kronecker product expression
//
//  Copyright (C) 2012-2020 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. You can redistribute it and/or modify it under
//  the terms of the New (Revised) BSD License. Redistribution and use in source and binary
//  forms, with or without modification, are permitted provided that the following conditions
//  are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice, this list
//     of conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//  3. Neither the names of the Blaze development group nor the names of its contributors
//     may be used to endorse or promote products derived from this software without specific
//     prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
//  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
//  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
//  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
//  DAMAGE.
*/
//=================================================================================================

#ifndef _BLAZE_MATH_EXPRESSIONS_TSMATSMATKRONEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_TSMATSMATKRONEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/Identity.h>
#include <blaze/math/constraints/MatMatKronExpr.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/constraints/Zero.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/MatMatKronExpr.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/traits/KronTrait.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsIdentity.h>
#include <blaze/math/typetraits/IsTemporary.h>
#include <blaze/math/typetraits/IsZero.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS TSMATSMATKRONEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for sparse matrix-sparse matrix Kronecker product.
// \ingroup sparse_matrix_expression
//
// The TSMatSMatKronExpr class represents the compile time expression for Kronecker products
// between a column-major sparse matrix and a row-major sparse matrix.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
class TSMatSMatKronExpr
   : public MatMatKronExpr< SparseMatrix< TSMatSMatKronExpr<MT1,MT2>, false > >
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   using RT1 = ResultType_t<MT1>;     //!< Result type of the left-hand side sparse matrix expression.
   using RT2 = ResultType_t<MT2>;     //!< Result type of the right-hand side sparse matrix expression.
   using RN1 = ReturnType_t<MT1>;     //!< Return type of the left-hand side sparse matrix expression.
   using RN2 = ReturnType_t<MT2>;     //!< Return type of the right-hand side sparse matrix expression.
   using CT1 = CompositeType_t<MT1>;  //!< Composite type of the left-hand side sparse matrix expression.
   using CT2 = CompositeType_t<MT2>;  //!< Composite type of the right-hand side sparse matrix expression.
   //**********************************************************************************************

   //**Return type evaluation**********************************************************************
   //! Compilation switch for the selection of the subscript operator return type.
   /*! The \a returnExpr compile time constant expression is a compilation switch for the
       selection of the \a ReturnType. If either matrix operand returns a temporary vector
       or matrix, \a returnExpr will be set to \a false and the subscript operator will
       return it's result by value. Otherwise \a returnExpr will be set to \a true and
       the subscript operator may return it's result as an expression. */
   static constexpr bool returnExpr = ( !IsTemporary_v<RN1> && !IsTemporary_v<RN2> );

   //! Expression return type for the subscript operator.
   using ExprReturnType = decltype( std::declval<RN1>() * std::declval<RN2>() );
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this TSMatSMatKronExpr instance.
   using This = TSMatSMatKronExpr<MT1,MT2>;

   //! Base type of this TSMatSMatKronExpr instance.
   using BaseType = MatMatKronExpr< SparseMatrix<This,false> >;

   using ResultType    = KronTrait_t<RT1,RT2>;         //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<ResultType>;   //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<ResultType>;    //!< Resulting element type.

   //! Return type for expression template evaluations.
   using ReturnType = const If_t< returnExpr, ExprReturnType, ElementType >;

   //! Data type for composite expression templates.
   using CompositeType = const ResultType;

   //! Composite type of the left-hand side sparse matrix expression.
   using LeftOperand = If_t< IsExpression_v<MT1>, const MT1, const MT1& >;

   //! Composite type of the right-hand side sparse matrix expression.
   using RightOperand = If_t< IsExpression_v<MT2>, const MT2, const MT2& >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = false;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the TSMatSMatKronExpr class.
   //
   // \param lhs The left-hand side sparse matrix operand of the Kronecker product expression.
   // \param rhs The right-hand side sparse matrix operand of the Kronecker product expression.
   */
   inline TSMatSMatKronExpr( const MT1& lhs, const MT2& rhs ) noexcept
      : lhs_( lhs )  // Left-hand side sparse matrix of the Kronecker product expression
      , rhs_( rhs )  // Right-hand side sparse matrix of the Kronecker product expression
   {}
   //**********************************************************************************************

   //**Access operator*****************************************************************************
   /*!\brief 2D-access to the matrix elements.
   //
   // \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
   // \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator()( size_t i, size_t j ) const {
      BLAZE_INTERNAL_ASSERT( i < rows()   , "Invalid row access index"    );
      BLAZE_INTERNAL_ASSERT( j < columns(), "Invalid column access index" );
      return lhs_( i/rhs_.rows(), j/rhs_.columns() ) * rhs_( i%rhs_.rows(), j%rhs_.columns() );
   }
   //**********************************************************************************************

   //**At function*********************************************************************************
   /*!\brief Checked access to the matrix elements.
   //
   // \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
   // \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   // \exception std::out_of_range Invalid matrix access index.
   */
   inline ReturnType at( size_t i, size_t j ) const {
      if( i >= rows() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
      }
      if( j >= columns() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
      }
      return (*this)(i,j);
   }
   //**********************************************************************************************

   //**Rows function*******************************************************************************
   /*!\brief Returns the current number of rows of the matrix.
   //
   // \return The number of rows of the matrix.
   */
   inline size_t rows() const noexcept {
      return lhs_.rows() * rhs_.rows();
   }
   //**********************************************************************************************

   //**Columns function****************************************************************************
   /*!\brief Returns the current number of columns of the matrix.
   //
   // \return The number of columns of the matrix.
   */
   inline size_t columns() const noexcept {
      return lhs_.columns() * rhs_.columns();
   }
   //**********************************************************************************************

   //**NonZeros function***************************************************************************
   /*!\brief Returns the number of non-zero elements in the sparse matrix.
   //
   // \return The number of non-zero elements in the sparse matrix.
   */
   inline size_t nonZeros() const {
      return lhs_.nonZeros() * rhs_.nonZeros();
   }
   //**********************************************************************************************

   //**NonZeros function***************************************************************************
   /*!\brief Returns the number of non-zero elements in the specified row.
   //
   // \param i The index of the row.
   // \return The number of non-zero elements of row \a i.
   */
   inline size_t nonZeros( size_t i ) const {
      MAYBE_UNUSED( i );
      return 0UL;
   }
   //**********************************************************************************************

   //**Left operand access*************************************************************************
   /*!\brief Returns the left-hand side sparse matrix operand.
   //
   // \return The left-hand side sparse matrix operand.
   */
   inline LeftOperand leftOperand() const noexcept {
      return lhs_;
   }
   //**********************************************************************************************

   //**Right operand access************************************************************************
   /*!\brief Returns the right-hand side sparse matrix operand.
   //
   // \return The right-hand side sparse matrix operand.
   */
   inline RightOperand rightOperand() const noexcept {
      return rhs_;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can alias with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case the expression can alias effect is detected, \a false otherwise.
   */
   template< typename T >
   inline bool canAlias( const T* alias ) const noexcept {
      return ( lhs_.canAlias( alias ) || rhs_.canAlias( alias ) );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression is aliased with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case an alias effect is detected, \a false otherwise.
   */
   template< typename T >
   inline bool isAliased( const T* alias ) const noexcept {
      return ( lhs_.isAliased( alias ) || rhs_.isAliased( alias ) );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  lhs_;  //!< Left-hand side sparse matrix of the Kronecker product expression.
   RightOperand rhs_;  //!< Right-hand side sparse matrix of the Kronecker product expression.
   //**********************************************************************************************

   //**Assignment to dense matrices****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a transpose sparse matrix-sparse matrix Kronecker product to a dense
   //        matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Kronecker product expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose sparse matrix-
   // sparse matrix Kronecker product expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO2 >   // Storage order of the target dense matrix
   friend inline void assign( DenseMatrix<MT,SO2>& lhs, const TSMatSMatKronExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( rhs.rows() == 0UL || rhs.columns() == 0UL ) {
         return;
      }

      CT1 A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side sparse matrix operand
      CT2 B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );

      const size_t M( B.rows()    );
      const size_t N( B.columns() );

      for( size_t j=0UL; j<A.columns(); ++j )
         for( auto aelem=A.begin(j); aelem!=A.end(j); ++aelem )
            for( size_t k=0UL; k<M; ++k )
               for( auto belem=B.begin(k); belem!=B.end(k); ++belem )
                  (*lhs)(aelem->index()*M+k,j*N+belem->index()) = aelem->value() * belem->value();
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to row-major sparse matrices*****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a transpose sparse matrix-sparse matrix Kronecker product to a row-major
   //        sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side Kronecker product expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose sparse matrix-
   // sparse matrix Kronecker product expression to a row-major sparse matrix.
   */
   template< typename MT >  // Type of the target sparse matrix
   friend inline void assign( SparseMatrix<MT,false>& lhs, const TSMatSMatKronExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( rhs.rows() == 0UL || rhs.columns() == 0UL ) {
         return;
      }

      CT1 A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side sparse matrix operand
      CT2 B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );

      const size_t M( B.rows()    );
      const size_t N( B.columns() );

      // Counting the number of elements per row in A
      std::vector<size_t> nonzeros( A.rows(), 0UL );
      for( size_t j=0UL; j<A.columns(); ++j ) {
         const auto end( A.end(j) );
         for( auto aelem=A.begin(j); aelem!=end; ++aelem ) {
            ++nonzeros[aelem->index()];
         }
      }

      // Resizing the left-hand side sparse matrix
      for( size_t i=0UL; i<A.rows(); ++i ) {
         for( size_t j=0UL; j<M; ++j ) {
            (*lhs).reserve( i*M+j, nonzeros[i]*B.nonZeros(j) );
         }
      }

      // Performing the Kronecker product
      for( size_t j=0UL; j<A.columns(); ++j )
         for( auto aelem=A.begin(j); aelem!=A.end(j); ++aelem )
            for( size_t k=0UL; k<M; ++k )
               for( auto belem=B.begin(k); belem!=B.end(k); ++belem )
                  (*lhs).append( aelem->index()*M+k, j*N+belem->index(), aelem->value() * belem->value(), true );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to column-major sparse matrices**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a transpose sparse matrix-sparse matrix Kronecker product to a
   //        column-major sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side Kronecker product expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose sparse matrix-
   // sparse matrix Kronecker product expression to a column-major sparse matrix.
   */
   template< typename MT >  // Type of the target sparse matrix
   friend inline void assign( SparseMatrix<MT,true>& lhs, const TSMatSMatKronExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( rhs.rows() == 0UL || rhs.columns() == 0UL ) {
         return;
      }

      CT1 A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side sparse matrix operand
      CT2 B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );

      const size_t M( B.rows()    );
      const size_t N( B.columns() );

      // Counting the number of elements per column in B
      std::vector<size_t> nonzeros( B.columns(), 0UL );
      for( size_t i=0UL; i<B.rows(); ++i ) {
         const auto end( B.end(i) );
         for( auto belem=B.begin(i); belem!=end; ++belem ) {
            ++nonzeros[belem->index()];
         }
      }

      // Resizing the left-hand side sparse matrix
      for( size_t i=0UL; i<A.columns(); ++i ) {
         for( size_t j=0UL; j<N; ++j ) {
            (*lhs).reserve( i*N+j, A.nonZeros(i)*nonzeros[j] );
         }
      }

      // Performing the Kronecker product
      for( size_t j=0UL; j<A.columns(); ++j )
         for( auto aelem=A.begin(j); aelem!=A.end(j); ++aelem )
            for( size_t k=0UL; k<M; ++k )
               for( auto belem=B.begin(k); belem!=B.end(k); ++belem )
                  (*lhs).append( aelem->index()*M+k, j*N+belem->index(), aelem->value() * belem->value(), true );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense matrices*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a transpose sparse matrix-sparse matrix Kronecker product to
   //        a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Kronecker product expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a transpose sparse
   // matrix-sparse matrix Kronecker product expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO2 >   // Storage order of the target dense matrix
   friend inline void addAssign( DenseMatrix<MT,SO2>& lhs, const TSMatSMatKronExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( rhs.rows() == 0UL || rhs.columns() == 0UL ) {
         return;
      }

      CT1 A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side sparse matrix operand
      CT2 B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );

      const size_t M( B.rows()    );
      const size_t N( B.columns() );

      for( size_t j=0UL; j<A.columns(); ++j )
         for( auto aelem=A.begin(j); aelem!=A.end(j); ++aelem )
            for( size_t k=0UL; k<M; ++k )
               for( auto belem=B.begin(k); belem!=B.end(k); ++belem )
                  (*lhs)(aelem->index()*M+k,j*N+belem->index()) += aelem->value() * belem->value();
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse matrices******************************************************
   // No special implementation for the addition assignment to sparse matrices.
   //**********************************************************************************************

   //**Subtraction assignment to dense matrices****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a transpose sparse matrix-sparse matrix Kronecker product
   //        to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Kronecker product expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a transpose
   // sparse matrix-sparse matrix Kronecker product expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO2 >   // Storage order of the target dense matrix
   friend inline void subAssign( DenseMatrix<MT,SO2>& lhs, const TSMatSMatKronExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( rhs.rows() == 0UL || rhs.columns() == 0UL ) {
         return;
      }

      CT1 A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side sparse matrix operand
      CT2 B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );

      const size_t M( B.rows()    );
      const size_t N( B.columns() );

      for( size_t j=0UL; j<A.columns(); ++j )
         for( auto aelem=A.begin(j); aelem!=A.end(j); ++aelem )
            for( size_t k=0UL; k<M; ++k )
               for( auto belem=B.begin(k); belem!=B.end(k); ++belem )
                  (*lhs)(aelem->index()*M+k,j*N+belem->index()) -= aelem->value() * belem->value();
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse matrices***************************************************
   // No special implementation for the subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**Schur product assignment to dense matrices**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Schur product assignment of a transpose sparse matrix-sparse matrix Kronecker product
   //        to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Kronecker product expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a transpose
   // sparse matrix-sparse matrix Kronecker product expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO2 >   // Storage order of the target dense matrix
   friend inline void schurAssign( DenseMatrix<MT,SO2>& lhs, const TSMatSMatKronExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( rhs.rows() == 0UL || rhs.columns() == 0UL ) {
         return;
      }

      CT1 A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side sparse matrix operand
      CT2 B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );

      const size_t M( B.rows()    );
      const size_t N( B.columns() );

      for( size_t j=0UL; j<A.columns(); ++j )
      {
         size_t i( 0UL );

         for( auto aelem=A.begin(j); aelem!=A.end(j); ++aelem, ++i )
         {
            for( ; i<aelem->index(); ++i ) {
               for( size_t k=0UL; k<M; ++k )
                  for( size_t l=0UL; l<N; ++l )
                     reset( (*lhs)(i*M+k,j*N+l) );
            }

            for( size_t k=0UL; k<M; ++k )
            {
               size_t l( 0UL );

               for( auto belem=B.begin(k); belem!=B.end(k); ++belem, ++l ) {
                  for( ; l<belem->index(); ++l )
                     reset( (*lhs)(i*M+k,j*N+l) );
                  (*lhs)(i*M+k,j*N+l) *= aelem->value() * belem->value();
               }

               for( ; l<N; ++l ) {
                  reset( (*lhs)(i*M+k,j*N+l) );
               }
            }
         }

         for( ; i<A.rows(); ++i ) {
            for( size_t k=0UL; k<M; ++k )
               for( size_t l=0UL; l<N; ++l )
                  reset( (*lhs)(i*M+k,j*N+l) );
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Schur product assignment to sparse matrices*************************************************
   // No special implementation for the Schur product assignment to sparse matrices.
   //**********************************************************************************************

   //**Multiplication assignment to dense matrices*************************************************
   // No special implementation for the multiplication assignment to dense matrices.
   //**********************************************************************************************

   //**Multiplication assignment to sparse matrices************************************************
   // No special implementation for the multiplication assignment to sparse matrices.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_FORM_VALID_MATMATKRONEXPR( MT1, MT2 );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the Kronecker product between two sparse matrices
//        (\f$ A=B \otimes C \f$).
// \ingroup sparse_matrix
//
// \param lhs The left-hand side sparse matrix for the Kronecker product.
// \param rhs The right-hand side sparse matrix for the Kronecker product.
// \return The Kronecker product of the two matrices.
//
// This function implements a performance optimized treatment of the Kronecker product between
// two sparse matrices.
*/
template< typename MT1  // Type of the left-hand side sparse matrix
        , typename MT2  // Type of the right-hand side sparse matrix
        , DisableIf_t< ( IsIdentity_v<MT1> && IsIdentity_v<MT2> ) ||
                       ( IsZero_v<MT1> || IsZero_v<MT2> ) >* = nullptr >
inline const TSMatSMatKronExpr<MT1,MT2>
   tsmatsmatkron( const SparseMatrix<MT1,true>& lhs, const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return TSMatSMatKronExpr<MT1,MT2>( *lhs, *rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the Kronecker product between a two identity matrices
//        (\f$ A=B \otimes C \f$).
// \ingroup sparse_matrix
//
// \param lhs The left-hand side sparse matrix for the Kronecker product.
// \param rhs The right-hand side sparse matrix for the Kronecker product.
// \return The Kronecker product of the two matrices.
//
// This function implements a performance optimized treatment of the Kronecker product between
// two identity matrices. It returns an identity matrix.
*/
template< typename MT1  // Type of the left-hand side sparse matrix
        , typename MT2  // Type of the right-hand side sparse matrix
        , EnableIf_t< IsIdentity_v<MT1> && IsIdentity_v<MT2> >* = nullptr >
inline decltype(auto)
   tsmatsmatkron( const SparseMatrix<MT1,true>& lhs, const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const KronTrait_t< ResultType_t<MT1>, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( ReturnType );
   BLAZE_CONSTRAINT_MUST_BE_IDENTITY_MATRIX_TYPE( ReturnType );

   return ReturnType( (*lhs).rows()*(*rhs).rows() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the Kronecker product between a zero matrix and sparse matrix
//        (\f$ A=B \otimes C \f$).
// \ingroup sparse_matrix
//
// \param lhs The left-hand side sparse matrix for the Kronecker product.
// \param rhs The right-hand side sparse matrix for the Kronecker product.
// \return The Kronecker product of the two matrices.
//
// This function implements a performance optimized treatment of the Kronecker product between a
// zero matrix and a sparse matrix. It returns a zero matrix.
*/
template< typename MT1  // Type of the left-hand side sparse matrix
        , typename MT2  // Type of the right-hand side sparse matrix
        , EnableIf_t< IsZero_v<MT1> || IsZero_v<MT2> >* = nullptr >
inline decltype(auto)
   tsmatsmatkron( const SparseMatrix<MT1,true>& lhs, const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const KronTrait_t< ResultType_t<MT1>, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_ZERO_TYPE( ReturnType );

   return ReturnType( (*lhs).rows()*(*rhs).rows(), (*lhs).columns()*(*rhs).columns() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Operator for the Kronecker product of a column-major and a row-major sparse matrix
//        (\f$ A=B \otimes C \f$).
// \ingroup sparse_matrix
//
// \param lhs The left-hand side sparse matrix for the Kronecker product.
// \param rhs The right-hand side sparse matrix for the Kronecker product.
// \return The Kronecker product of the two matrices.
//
// The kron() function computes the Kronecker product for a column-major sparse matrix and a
// row-major sparse matrix:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   blaze::CompressedMatrix<double,columnMajor> A;
   blaze::CompressedMatrix<double,rowMajor> B, C;
   // ... Resizing and initialization
   C = kron( A, B );
   \endcode

// The function returns an expression representing a sparse matrix of the higher-order element
// type of the two involved matrix element types \a MT1::ElementType and \a MT2::ElementType.
// Both matrix types \a MT1 and \a MT2 as well as the two element types \a MT1::ElementType
// and \a MT2::ElementType have to be supported by the MultTrait class template.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
inline decltype(auto)
   kron( const SparseMatrix<MT1,true>& lhs, const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return tsmatsmatkron( *lhs, *rhs );
}
//*************************************************************************************************

} // namespace blaze

#endif
