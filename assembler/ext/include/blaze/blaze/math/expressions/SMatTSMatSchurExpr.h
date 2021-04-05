//=================================================================================================
/*!
//  \file blaze/math/expressions/SMatTSMatSchurExpr.h
//  \brief Header file for the sparse matrix/transpose sparse matrix Schur product expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_SMATTSMATSCHUREXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_SMATTSMATSCHUREXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/Identity.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
#include <blaze/math/constraints/SchurExpr.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/Symmetric.h>
#include <blaze/math/constraints/Zero.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/SchurExpr.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/sparse/Forward.h>
#include <blaze/math/traits/SchurTrait.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/typetraits/IsTemporary.h>
#include <blaze/math/typetraits/IsUniLower.h>
#include <blaze/math/typetraits/IsUniUpper.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/typetraits/IsZero.h>
#include <blaze/math/typetraits/StorageOrder.h>
#include <blaze/util/algorithms/Min.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS SMATTSMATSCHUREXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for sparse matrix-transpose sparse matrix Schur product.
// \ingroup sparse_matrix_expression
//
// The SMatTSMatSchurExpr class represents the compile time expression for Schur products between
// a row-major sparse matrix and a column-major sparse matrix.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
class SMatTSMatSchurExpr
   : public SchurExpr< SparseMatrix< SMatTSMatSchurExpr<MT1,MT2>, false > >
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

   //**Serial evaluation strategy******************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! This variable template is a helper for the selection of the serial evaluation strategy.
       In case the two given matrix types have a different storage order and in case the second
       matrix type is symmetric, the variable is set to 1 and an optimized evaluation strategy
       is selected. Otherwise the variable is set to 0 and the default strategy is chosen. */
   template< typename T1, typename T2 >
   static constexpr bool UseSymmetricKernel_v =
      ( ( StorageOrder_v<T1> != StorageOrder_v<T2> ) && IsSymmetric_v<T2> );
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this SMatTSMatSchurExpr instance.
   using This = SMatTSMatSchurExpr<MT1,MT2>;

   //! Base type of this SMatTSMatSchurExpr instance.
   using BaseType = SchurExpr< SparseMatrix<This,false> >;

   using ResultType    = SchurTrait_t<RT1,RT2>;        //!< Result type for expression template evaluations.
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
   /*!\brief Constructor for the SMatTSMatSchurExpr class.
   //
   // \param lhs The left-hand side sparse matrix operand of the Schur product expression.
   // \param rhs The right-hand side sparse matrix operand of the Schur product expression.
   */
   inline SMatTSMatSchurExpr( const MT1& lhs, const MT2& rhs ) noexcept
      : lhs_( lhs )  // Left-hand side sparse matrix of the Schur product expression
      , rhs_( rhs )  // Right-hand side sparse matrix of the Schur product expression
   {
      BLAZE_INTERNAL_ASSERT( lhs.rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( lhs.columns() == rhs.columns(), "Invalid number of columns" );
   }
   //**********************************************************************************************

   //**Access operator*****************************************************************************
   /*!\brief 2D-access to the matrix elements.
   //
   // \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
   // \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator()( size_t i, size_t j ) const {
      BLAZE_INTERNAL_ASSERT( i < lhs_.rows()   , "Invalid row access index"    );
      BLAZE_INTERNAL_ASSERT( j < lhs_.columns(), "Invalid column access index" );
      return lhs_(i,j) * rhs_(i,j);
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
      if( i >= lhs_.rows() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
      }
      if( j >= lhs_.columns() ) {
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
      return lhs_.rows();
   }
   //**********************************************************************************************

   //**Columns function****************************************************************************
   /*!\brief Returns the current number of columns of the matrix.
   //
   // \return The number of columns of the matrix.
   */
   inline size_t columns() const noexcept {
      return lhs_.columns();
   }
   //**********************************************************************************************

   //**NonZeros function***************************************************************************
   /*!\brief Returns the number of non-zero elements in the sparse matrix.
   //
   // \return The number of non-zero elements in the sparse matrix.
   */
   inline size_t nonZeros() const {
      return min( lhs_.nonZeros(), rhs_.nonZeros() );
   }
   //**********************************************************************************************

   //**NonZeros function***************************************************************************
   /*!\brief Returns the number of non-zero elements in the specified row.
   //
   // \param i The index of the row.
   // \return The number of non-zero elements of row \a i.
   */
   inline size_t nonZeros( size_t i ) const {
      return min( lhs_.nonZeros(i), rhs_.nonZeros(i) );
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
   /*!\brief Returns the right-hand side transpose sparse matrix operand.
   //
   // \return The right-hand side transpose sparse matrix operand.
   */
   inline RightOperand rightOperand() const noexcept {
      return rhs_;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can alias with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case the expression can alias, \a false otherwise.
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
   LeftOperand  lhs_;  //!< Left-hand side sparse matrix of the Schur product expression.
   RightOperand rhs_;  //!< Right-hand side sparse matrix of the Schur product expression.
   //**********************************************************************************************

   //**Assignment to row-major dense matrices******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse matrix-transpose sparse matrix Schur product to a row-major
   //        dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Schur product expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse matrix-transpose
   // sparse matrix Schur product expression to a row-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline auto assign( DenseMatrix<MT,false>& lhs, const SMatTSMatSchurExpr& rhs )
      -> DisableIf_t< UseSymmetricKernel_v<MT,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      // Evaluation of the left-hand side sparse matrix operand
      CT1 A( serial( rhs.lhs_ ) );

      // Evaluation of the right-hand side sparse matrix operand
      const OppositeType_t<RT2> B( serial( rhs.rhs_ ) );

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == (*lhs).columns()  , "Invalid number of columns" );

      assign( *lhs, A % B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to column-major dense matrices***************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse matrix-transpose sparse matrix Schur product to a column-major
   //        dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Schur product expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse matrix-transpose
   // sparse matrix Schur product expression to a column-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline auto assign( DenseMatrix<MT,true>& lhs, const SMatTSMatSchurExpr& rhs )
      -> DisableIf_t< UseSymmetricKernel_v<MT,MT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      // Evaluation of the left-hand side sparse matrix operand
      const OppositeType_t<RT1> A( serial( rhs.lhs_ ) );

      // Evaluation of the right-hand side sparse matrix operand
      CT2 B( serial( rhs.rhs_ ) );

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == (*lhs).columns()  , "Invalid number of columns" );

      assign( *lhs, A % B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to row-major sparse matrices*****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse matrix-transpose sparse matrix Schur product to a row-major
   //        sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side Schur product expression to be assigned.
   // \return void
   //
   // This function implements the default assignment of a sparse matrix-transpose sparse matrix
   // Schur product expression to a row-major sparse matrix.
   */
   template< typename MT >  // Type of the target sparse matrix
   friend inline auto assign( SparseMatrix<MT,false>& lhs, const SMatTSMatSchurExpr& rhs )
      -> DisableIf_t< UseSymmetricKernel_v<MT,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      // Evaluation of the left-hand side sparse matrix operand
      CT1 A( serial( rhs.lhs_ ) );

      // Evaluation of the right-hand side sparse matrix operand
      const OppositeType_t<RT2> B( serial( rhs.rhs_ ) );

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == (*lhs).columns()  , "Invalid number of columns" );

      assign( *lhs, A % B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to column-major sparse matrices**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse matrix-transpose sparse matrix Schur product to a column-major
   //        sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side Schur product expression to be assigned.
   // \return void
   //
   // This function implements the default assignment of a sparse matrix-transpose sparse matrix
   // Schur product expression to a column-major sparse matrix.
   */
   template< typename MT >  // Type of the target sparse matrix
   friend inline auto assign( SparseMatrix<MT,true>& lhs, const SMatTSMatSchurExpr& rhs )
      -> DisableIf_t< UseSymmetricKernel_v<MT,MT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      // Evaluation of the left-hand side sparse matrix operand
      const OppositeType_t<RT1> A( serial( rhs.lhs_ ) );

      // Evaluation of the right-hand side sparse matrix operand
      CT2 B( serial( rhs.rhs_ ) );

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == (*lhs).columns()  , "Invalid number of columns" );

      assign( *lhs, A % B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Restructuring assignment to row-major matrices**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Restructuring assignment of a sparse matrix-transpose sparse matrix Schur product to
   //        a row-major matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side Schur product expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse matrix-transpose
   // sparse matrix Schur product expression to a row-major matrix.
   */
   template< typename MT >  // Type of the target sparse matrix
   friend inline auto assign( Matrix<MT,false>& lhs, const SMatTSMatSchurExpr& rhs )
      -> EnableIf_t< UseSymmetricKernel_v<MT,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      assign( *lhs, rhs.lhs_ % trans( rhs.rhs_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Restructuring assignment to column-major matrices*******************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Restructuring assignment of a sparse matrix-transpose sparse matrix Schur product to
   //        a column-major sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side Schur product expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse matrix-transpose
   // sparse matrix Schur product expression to a column-major matrix.
   */
   template< typename MT >  // Type of the target sparse matrix
   friend inline auto assign( Matrix<MT,true>& lhs, const SMatTSMatSchurExpr& rhs )
      -> EnableIf_t< UseSymmetricKernel_v<MT,MT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      assign( *lhs, trans( rhs.lhs_ ) % rhs.rhs_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to row-major dense matrices*********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a sparse matrix-transpose sparse matrix Schur product to a
   //        row-major dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Schur product expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a sparse
   // matrix-transpose sparse matrix Schur product expression to a row-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline auto addAssign( DenseMatrix<MT,false>& lhs, const SMatTSMatSchurExpr& rhs )
      -> DisableIf_t< UseSymmetricKernel_v<MT,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      // Evaluation of the left-hand side sparse matrix operand
      CT1 A( serial( rhs.lhs_ ) );

      // Evaluation of the right-hand side sparse matrix operand
      const OppositeType_t<RT2> B( serial( rhs.rhs_ ) );

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == (*lhs).columns()  , "Invalid number of columns" );

      addAssign( *lhs, A % B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to column-major dense matrices******************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a sparse matrix-transpose sparse matrix Schur product to a
   //        column-major dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Schur product expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a sparse
   // matrix-transpose sparse matrix Schur product expression to a column-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline auto addAssign( DenseMatrix<MT,true>& lhs, const SMatTSMatSchurExpr& rhs )
      -> DisableIf_t< UseSymmetricKernel_v<MT,MT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      // Evaluation of the left-hand side sparse matrix operand
      const OppositeType_t<RT1> A( serial( rhs.lhs_ ) );

      // Evaluation of the right-hand side sparse matrix operand
      CT2 B( serial( rhs.rhs_ ) );

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == (*lhs).columns()  , "Invalid number of columns" );

      addAssign( *lhs, A % B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Restructuring addition assignment to row-major matrices*************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Restructuring addition assignment of a sparse matrix-transpose sparse matrix Schur
   //        product to a row-major matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side Schur product expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a sparse matrix-
   // transpose sparse matrix Schur product expression to a row-major matrix.
   */
   template< typename MT >  // Type of the target sparse matrix
   friend inline auto addAssign( Matrix<MT,false>& lhs, const SMatTSMatSchurExpr& rhs )
      -> EnableIf_t< UseSymmetricKernel_v<MT,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      addAssign( *lhs, rhs.lhs_ % trans( rhs.rhs_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Restructuring addition assignment to column-major matrices**********************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Restructuring addition assignment of a sparse matrix-transpose sparse matrix Schur
   //        product to a column-major sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side Schur product expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a sparse matrix-
   // transpose sparse matrix Schur product expression to a column-major matrix.
   */
   template< typename MT >  // Type of the target sparse matrix
   friend inline auto addAssign( Matrix<MT,true>& lhs, const SMatTSMatSchurExpr& rhs )
      -> EnableIf_t< UseSymmetricKernel_v<MT,MT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      addAssign( *lhs, trans( rhs.lhs_ ) % rhs.rhs_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse matrices******************************************************
   // No special implementation for the addition assignment to sparse matrices.
   //**********************************************************************************************

   //**Subtraction assignment to row-major dense matrices******************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a sparse matrix-transpose sparse matrix Schur product to a
   //        row-major dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Schur product expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a sparse
   // matrix-transpose sparse matrix Schur product expression to a row-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline auto subAssign( DenseMatrix<MT,false>& lhs, const SMatTSMatSchurExpr& rhs )
      -> DisableIf_t< UseSymmetricKernel_v<MT,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      // Evaluation of the left-hand side sparse matrix operand
      CT1 A( serial( rhs.lhs_ ) );

      // Evaluation of the right-hand side sparse matrix operand
      const OppositeType_t<RT2> B( serial( rhs.rhs_ ) );

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == (*lhs).columns()  , "Invalid number of columns" );

      subAssign( *lhs, A % B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to column-major dense matrices***************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a sparse matrix-transpose sparse matrix Schur product to a
   //        column-major dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Schur product expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a sparse
   // matrix-transpose sparse matrix Schur product expression to a column-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline auto subAssign( DenseMatrix<MT,true>& lhs, const SMatTSMatSchurExpr& rhs )
      -> DisableIf_t< UseSymmetricKernel_v<MT,MT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      // Evaluation of the left-hand side sparse matrix operand
      const OppositeType_t<RT1> A( serial( rhs.lhs_ ) );

      // Evaluation of the right-hand side sparse matrix operand
      CT2 B( serial( rhs.rhs_ ) );

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == (*lhs).columns()  , "Invalid number of columns" );

      subAssign( *lhs, A % B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Restructuring subtraction assignment to row-major matrices**********************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Restructuring subtraction assignment of a sparse matrix-transpose sparse matrix Schur
   //        product to a row-major matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side Schur product expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a sparse matrix-
   // transpose sparse matrix Schur product expression to a row-major matrix.
   */
   template< typename MT >  // Type of the target sparse matrix
   friend inline auto subAssign( Matrix<MT,false>& lhs, const SMatTSMatSchurExpr& rhs )
      -> EnableIf_t< UseSymmetricKernel_v<MT,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      subAssign( *lhs, rhs.lhs_ % trans( rhs.rhs_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Restructuring subtraction assignment to column-major matrices*******************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Restructuring subtraction assignment of a sparse matrix-transpose sparse matrix Schur
   //        product to a column-major sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side Schur product expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a sparse matrix-
   // transpose sparse matrix Schur product expression to a column-major matrix.
   */
   template< typename MT >  // Type of the target sparse matrix
   friend inline auto subAssign( Matrix<MT,true>& lhs, const SMatTSMatSchurExpr& rhs )
      -> EnableIf_t< UseSymmetricKernel_v<MT,MT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      subAssign( *lhs, trans( rhs.lhs_ ) % rhs.rhs_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse matrices***************************************************
   // No special implementation for the subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**Schur product assignment to row-major dense matrices****************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Schur product assignment of a sparse matrix-transpose sparse matrix Schur product to
   //        a row-major dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Schur product expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a sparse
   // matrix-transpose sparse matrix Schur product expression to a row-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline auto schurAssign( DenseMatrix<MT,false>& lhs, const SMatTSMatSchurExpr& rhs )
      -> DisableIf_t< UseSymmetricKernel_v<MT,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      // Evaluation of the left-hand side sparse matrix operand
      CT1 A( serial( rhs.lhs_ ) );

      // Evaluation of the right-hand side sparse matrix operand
      const OppositeType_t<RT2> B( serial( rhs.rhs_ ) );

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == (*lhs).columns()  , "Invalid number of columns" );

      schurAssign( *lhs, A % B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Schur product assignment to column-major dense matrices*************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Schur product assignment of a sparse matrix-transpose sparse matrix Schur product to
   //        a column-major dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Schur product expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a sparse
   // matrix-transpose sparse matrix Schur product expression to a column-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline auto schurAssign( DenseMatrix<MT,true>& lhs, const SMatTSMatSchurExpr& rhs )
      -> DisableIf_t< UseSymmetricKernel_v<MT,MT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      // Evaluation of the left-hand side sparse matrix operand
      const OppositeType_t<RT1> A( serial( rhs.lhs_ ) );

      // Evaluation of the right-hand side sparse matrix operand
      CT2 B( serial( rhs.rhs_ ) );

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == (*lhs).columns()  , "Invalid number of columns" );

      schurAssign( *lhs, A % B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Restructuring Schur product assignment to row-major matrices********************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Restructuring Schur product assignment of a sparse matrix-transpose sparse matrix
   //        Schur product to a row-major matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side Schur product expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a sparse
   // matrix-transpose sparse matrix Schur product expression to a row-major matrix.
   */
   template< typename MT >  // Type of the target sparse matrix
   friend inline auto schurAssign( Matrix<MT,false>& lhs, const SMatTSMatSchurExpr& rhs )
      -> EnableIf_t< UseSymmetricKernel_v<MT,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      schurAssign( *lhs, rhs.lhs_ % trans( rhs.rhs_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Restructuring Schur product assignment to column-major matrices*****************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Restructuring Schur product assignment of a sparse matrix-transpose sparse matrix
   //        Schur product to a column-major sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side Schur product expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a sparse
   // matrix-transpose sparse matrix Schur product expression to a column-major matrix.
   */
   template< typename MT >  // Type of the target sparse matrix
   friend inline auto schurAssign( Matrix<MT,true>& lhs, const SMatTSMatSchurExpr& rhs )
      -> EnableIf_t< UseSymmetricKernel_v<MT,MT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      schurAssign( *lhs, trans( rhs.lhs_ ) % rhs.rhs_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Multiplication assignment to dense matrices*************************************************
   // No special implementation for the multiplication assignment to dense matrices.
   //**********************************************************************************************

   //**Multiplication assignment to sparse matrices************************************************
   // No special implementation for the multiplication assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP assignment to dense matrices************************************************************
   // No special implementation for the SMP assignment to dense matrices.
   //**********************************************************************************************

   //**SMP assignment to sparse matrices***********************************************************
   // No special implementation for the SMP assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP addition assignment to dense matrices***************************************************
   // No special implementation for the SMP addition assignment to dense matrices.
   //**********************************************************************************************

   //**SMP addition assignment to sparse matrices**************************************************
   // No special implementation for the SMP addition assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP subtraction assignment to dense matrices************************************************
   // No special implementation for the SMP subtraction assignment to dense matrices.
   //**********************************************************************************************

   //**SMP subtraction assignment to sparse matrices***********************************************
   // No special implementation for the SMP subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP Schur product assignment to dense matrices**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP Schur product assignment of a sparse matrix-transpose sparse matrix Schur product
   //        to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Schur product expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized SMP Schur product assignment of a
   // sparse matrix-transpose sparse matrix Schur product expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline void smpSchurAssign( DenseMatrix<MT,SO>& lhs, const SMatTSMatSchurExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      smpSchurAssign( *lhs, rhs.lhs_ );
      smpSchurAssign( *lhs, rhs.rhs_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP Schur product assignment to sparse matrices*********************************************
   // No special implementation for the SMP Schur product assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP multiplication assignment to dense matrices*********************************************
   // No special implementation for the SMP multiplication assignment to dense matrices.
   //**********************************************************************************************

   //**SMP multiplication assignment to sparse matrices********************************************
   // No special implementation for the SMP multiplication assignment to sparse matrices.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_FORM_VALID_SCHUREXPR( MT1, MT2 );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the Schur product between a row-major sparse matrix and a
//        column-major sparse matrix (\f$ A=B \circ C \f$).
// \ingroup sparse_matrix
//
// \param lhs The left-hand side sparse matrix for the Schur product.
// \param rhs The right-hand side sparse matrix for the Schur product.
// \return The Schur product of the two matrices.
//
// This function implements a performance optimized treatment of the Schur product between a
// row-major sparse matrix and a column-major sparse matrix.
*/
template< typename MT1  // Type of the left-hand side sparse matrix
        , typename MT2  // Type of the right-hand side sparse matrix
        , DisableIf_t< ( IsUniLower_v<MT1> && IsUniUpper_v<MT2> ) ||
                       ( IsUniUpper_v<MT1> && IsUniLower_v<MT2> ) ||
                       ( IsStrictlyLower_v<MT1> && IsUpper_v<MT2> ) ||
                       ( IsStrictlyUpper_v<MT1> && IsLower_v<MT2> ) ||
                       ( IsLower_v<MT1> && IsStrictlyUpper_v<MT2> ) ||
                       ( IsUpper_v<MT1> && IsStrictlyLower_v<MT2> ) ||
                       ( IsZero_v<MT1> || IsZero_v<MT2> ) >* = nullptr >
inline const SMatTSMatSchurExpr<MT1,MT2>
   smattsmatschur( const SparseMatrix<MT1,false>& lhs, const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).columns(), "Invalid number of columns" );

   return SMatTSMatSchurExpr<MT1,MT2>( *lhs, *rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the Schur product between a unitriangular row-major sparse
//        matrix and a unitriangular column-major sparse matrix (\f$ A=B \circ C \f$).
// \ingroup sparse_matrix
//
// \param lhs The left-hand side sparse matrix for the Schur product.
// \param rhs The right-hand side sparse matrix for the Schur product.
// \return The resulting identity matrix.
//
// This function implements a performance optimized treatment of the Schur product between a
// unitriangular row-major sparse matrix and a unitriangular column-major sparse matrix. It
// returns an identity matrix.
*/
template< typename MT1  // Type of the left-hand side sparse matrix
        , typename MT2  // Type of the right-hand side sparse matrix
        , EnableIf_t< ( IsUniLower_v<MT1> && IsUniUpper_v<MT2> ) ||
                      ( IsUniUpper_v<MT1> && IsUniLower_v<MT2> ) >* = nullptr >
inline decltype(auto)
   smattsmatschur( const SparseMatrix<MT1,false>& lhs, const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( rhs );

   BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).columns(), "Invalid number of columns" );

   using ReturnType = const SchurTrait_t< ResultType_t<MT1>, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( ReturnType );
   BLAZE_CONSTRAINT_MUST_BE_IDENTITY_MATRIX_TYPE( ReturnType );

   return ReturnType( (*lhs).rows() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the Schur product between a (strictly) triangular row-major
//        sparse matrix and a (strictly) triangular column-major sparse matrix (\f$ A=B \circ C \f$).
// \ingroup sparse_matrix
//
// \param lhs The left-hand side sparse matrix for the Schur product.
// \param rhs The right-hand side sparse matrix for the Schur product.
// \return The resulting zero matrix.
//
// This function implements a performance optimized treatment of the Schur product between a
// (strictly) triangular row-major sparse matrix and a (strictly) triangular column-major sparse
// matrix. It returns a zero matrix.
*/
template< typename MT1  // Type of the left-hand side sparse matrix
        , typename MT2  // Type of the right-hand side sparse matrix
        , EnableIf_t< ( IsStrictlyLower_v<MT1> && IsUpper_v<MT2> ) ||
                      ( IsStrictlyUpper_v<MT1> && IsLower_v<MT2> ) ||
                      ( IsLower_v<MT1> && IsStrictlyUpper_v<MT2> ) ||
                      ( IsUpper_v<MT1> && IsStrictlyLower_v<MT2> ) ||
                      ( IsZero_v<MT1> || IsZero_v<MT2> ) >* = nullptr >
inline decltype(auto)
   smattsmatschur( const SparseMatrix<MT1,false>& lhs, const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( rhs );

   BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).columns(), "Invalid number of columns" );

   using ReturnType = const SchurTrait_t< ResultType_t<MT1>, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( ReturnType );
   BLAZE_CONSTRAINT_MUST_BE_ZERO_TYPE( ReturnType );

   return ReturnType( (*lhs).rows(), (*lhs).columns() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Operator for the Schur product of a row-major and a column-major sparse matrix
//        (\f$ A=B \circ C \f$).
// \ingroup sparse_matrix
//
// \param lhs The left-hand side sparse matrix for the Schur product.
// \param rhs The right-hand side sparse matrix for the Schur product.
// \return The Schur product of the two matrices.
// \exception std::invalid_argument Matrix sizes do not match.
//
// This operator represents the Schur product of a row-major and a column-major sparse matrix:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   blaze::CompressedMatrix<double,rowMajor> A, C;
   blaze::CompressedMatrix<double,columnMajor> B;
   // ... Resizing and initialization
   C = A % B;
   \endcode

// The operator returns an expression representing a sparse matrix of the higher-order element
// type of the two involved matrix element types \a MT1::ElementType and \a MT2::ElementType.
// Both matrix types \a MT1 and \a MT2 as well as the two element types \a MT1::ElementType
// and \a MT2::ElementType have to be supported by the MultTrait class template.\n
// In case the current number of rows and columns of the two given matrices don't match, a
// \a std::invalid_argument is thrown.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
inline decltype(auto)
   operator%( const SparseMatrix<MT1,false>& lhs, const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   if( (*lhs).rows() != (*rhs).rows() || (*lhs).columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   return smattsmatschur( *lhs, *rhs );
}
//*************************************************************************************************

} // namespace blaze

#endif
