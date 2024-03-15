//=================================================================================================
/*!
//  \file blaze/math/expressions/DMatTSMatAddExpr.h
//  \brief Header file for the dense matrix/transpose sparse matrix addition expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DMATTSMATADDEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DMATTSMATADDEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/MatMatAddExpr.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/MatMatAddExpr.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsOperation.h>
#include <blaze/math/typetraits/IsTemporary.h>
#include <blaze/math/typetraits/IsZero.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsSame.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DMATTSMATADDEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for dense matrix-sparse matrix additions.
// \ingroup dense_matrix_expression
//
// The DMatTSMatAddExpr class represents the compile time expression for additions between
// a row-major dense matrix and a column-major sparse matrix.
*/
template< typename MT1   // Type of the left-hand side dense matrix
        , typename MT2 > // Type of the right-hand side sparse matrix
class DMatTSMatAddExpr
   : public MatMatAddExpr< DenseMatrix< DMatTSMatAddExpr<MT1,MT2>, false > >
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
  using RT1 = ResultType_t<MT1>;  //!< Result type of the left-hand side dense matrix expression.
  using RT2 = ResultType_t<MT2>;  //!< Result type of the right-hand side sparse matrix expression.
  using RN1 = ReturnType_t<MT1>;  //!< Return type of the left-hand side dense matrix expression.
  using RN2 = ReturnType_t<MT2>;  //!< Return type of the right-hand side sparse matrix expression.
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
   using ExprReturnType = decltype( std::declval<RN1>() + std::declval<RN2>() );
   //**********************************************************************************************

   //**Parallel evaluation strategy****************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! This variable template is a helper for the selection of the parallel evaluation strategy.
       In case at least one of the two matrix operands is not SMP assignable, the variable is set
       to 1 and the expression specific evaluation strategy is selected. Otherwise the variable
       is set to 0 and the default strategy is chosen. */
   template< typename MT >
   static constexpr bool UseSMPAssign_v = ( !MT1::smpAssignable || !MT2::smpAssignable );
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this DMatTSMatAddExpr instance.
   using This = DMatTSMatAddExpr<MT1,MT2>;

   //! Base type of this DMatTSMatAddExpr instance.
   using BaseType = MatMatAddExpr< DenseMatrix<This,false> >;

   using ResultType    = AddTrait_t<RT1,RT2>;          //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<ResultType>;   //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<ResultType>;    //!< Resulting element type.

   //! Return type for expression template evaluations.
   using ReturnType = const If_t< returnExpr, ExprReturnType, ElementType >;

   //! Data type for composite expression templates.
   using CompositeType = const ResultType;

   //! Composite type of the left-hand side dense matrix expression.
   using LeftOperand = If_t< IsExpression_v<MT1>, const MT1, const MT1& >;

   //! Composite type of the right-hand side sparse matrix expression.
   using RightOperand = If_t< IsExpression_v<MT2>, const MT2, const MT2& >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = false;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = false;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DMatTSMatAddExpr class.
   //
   // \param lhs The left-hand side dense matrix operand of the addition expression.
   // \param rhs The right-hand side sparse matrix operand of the addition expression.
   */
   inline DMatTSMatAddExpr( const MT1& lhs, const MT2& rhs ) noexcept
      : lhs_( lhs )  // Left-hand side dense matrix of the addition expression
      , rhs_( rhs )  // Right-hand side sparse matrix of the addition expression
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
      return lhs_(i,j) + rhs_(i,j);
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

   //**Left operand access*************************************************************************
   /*!\brief Returns the left-hand side dense matrix operand.
   //
   // \return The left-hand side dense matrix operand.
   */
   inline LeftOperand leftOperand() const {
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
      return ( IsExpression_v<MT1> && lhs_.canAlias( alias ) ) ||
             ( rhs_.canAlias( alias ) );
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
      return ( lhs_.canAlias( alias ) || rhs_.canAlias( alias ) );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  lhs_;  //!< Left-hand side dense matrix of the addition expression.
   RightOperand rhs_;  //!< Right-hand side sparse matrix of the addition expression.
   //**********************************************************************************************

   //**Assignment to dense matrices****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense matrix-transpose sparse matrix addition to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side addition expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix-transpose
   // sparse matrix addition expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO2 >   // Storage order of the target dense matrix
   friend inline void assign( DenseMatrix<MT,SO2>& lhs, const DMatTSMatAddExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( !IsOperation_v<MT1> && isSame( *lhs, rhs.lhs_ ) ) {
         addAssign( *lhs, rhs.rhs_ );
      }
      else {
         assign   ( *lhs, rhs.lhs_ );
         addAssign( *lhs, rhs.rhs_ );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse matrices***************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense matrix-transpose sparse matrix addition to a sparse matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side addition expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix-transpose
   // sparse matrix addition expression to a sparse matrix.
   */
   template< typename MT  // Type of the target sparse matrix
           , bool SO2 >   // Storage order of the target sparse matrix
   friend inline void assign( SparseMatrix<MT,SO2>& lhs, const DMatTSMatAddExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      using TmpType = If_t< SO2, OppositeType, ResultType >;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( OppositeType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( OppositeType );
      BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER( MT, TmpType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( TmpType );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const TmpType tmp( serial( rhs ) );
      assign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense matrices*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a dense matrix-transpose sparse matrix addition to a dense
   //        matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side addition expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense matrix-
   // transpose sparse matrix addition expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO2 >   // Storage order of the target dense matrix
   friend inline void addAssign( DenseMatrix<MT,SO2>& lhs, const DMatTSMatAddExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      addAssign( *lhs, rhs.lhs_ );
      addAssign( *lhs, rhs.rhs_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse matrices******************************************************
   // No special implementation for the addition assignment to sparse matrices.
   //**********************************************************************************************

   //**Subtraction assignment to dense matrices****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a dense matrix-transpose sparse matrix addition to a dense
   //        matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side addition expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a dense matrix-
   // transpose sparse matrix addition expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO2 >   // Storage order of the target dense matrix
   friend inline void subAssign( DenseMatrix<MT,SO2>& lhs, const DMatTSMatAddExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      subAssign( *lhs, rhs.lhs_ );
      subAssign( *lhs, rhs.rhs_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse matrices***************************************************
   // No special implementation for the subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**Schur product assignment to dense matrices**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Schur product assignment of a dense matrix-transpose sparse matrix addition to a
   //        dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side addition expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a dense
   // matrix-transpose sparse matrix addition expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO2 >   // Storage order of the target dense matrix
   friend inline void schurAssign( DenseMatrix<MT,SO2>& lhs, const DMatTSMatAddExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ResultType tmp( serial( rhs ) );
      schurAssign( *lhs, tmp );
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

   //**SMP assignment to dense matrices************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a dense matrix-transpose sparse matrix addition to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side addition expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a dense matrix-transpose
   // sparse matrix addition expression to a dense matrix. Due to the explicit application of the
   // SFINAE principle, this function can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO2 >   // Storage order of the target dense matrix
   friend inline auto smpAssign( DenseMatrix<MT,SO2>& lhs, const DMatTSMatAddExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( !IsOperation_v<MT1> && isSame( *lhs, rhs.lhs_ ) ) {
         smpAddAssign( *lhs, rhs.rhs_ );
      }
      else {
         smpAssign   ( *lhs, rhs.lhs_ );
         smpAddAssign( *lhs, rhs.rhs_ );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to sparse matrices***********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a dense matrix-transpose sparse matrix addition to a sparse matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side addition expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a dense matrix-transpose
   // sparse matrix addition expression to a sparse matrix. Due to the explicit application of the
   // SFINAE principle, this function can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename MT  // Type of the target sparse matrix
           , bool SO2 >   // Storage order of the target sparse matrix
   friend inline auto smpAssign( SparseMatrix<MT,SO2>& lhs, const DMatTSMatAddExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT> >
   {
      BLAZE_FUNCTION_TRACE;

      using TmpType = If_t< SO2, OppositeType, ResultType >;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( OppositeType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( OppositeType );
      BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER( MT, TmpType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( TmpType );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const TmpType tmp( rhs );
      smpAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to dense matrices***************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a dense matrix-transpose sparse matrix addition to a dense
   //        matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side addition expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense matrix-
   // transpose sparse matrix addition expression to a dense matrix. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO2 >   // Storage order of the target dense matrix
   friend inline auto smpAddAssign( DenseMatrix<MT,SO2>& lhs, const DMatTSMatAddExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      smpAddAssign( *lhs, rhs.lhs_ );
      smpAddAssign( *lhs, rhs.rhs_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to sparse matrices**************************************************
   // No special implementation for the SMP addition assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP subtraction assignment to dense matrices************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a dense matrix-transpose sparse matrix addition to a dense
   //        matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side addition expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a dense
   // matrix-transpose sparse matrix addition expression to a dense matrix. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler in
   // case the expression specific parallel evaluation strategy is selected.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO2 >   // Storage order of the target dense matrix
   friend inline auto smpSubAssign( DenseMatrix<MT,SO2>& lhs, const DMatTSMatAddExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      smpSubAssign( *lhs, rhs.lhs_ );
      smpSubAssign( *lhs, rhs.rhs_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to sparse matrices***********************************************
   // No special implementation for the SMP subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP Schur product assignment to dense matrices**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP Schur product assignment of a dense matrix-transpose sparse matrix addition to
   //        a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side addition expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized SMP Schur product assignment of a dense
   // matrix-transpose sparse matrix addition expression to a dense matrix. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler in
   // case the expression specific parallel evaluation strategy is selected.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO2 >   // Storage order of the target dense matrix
   friend inline auto smpSchurAssign( DenseMatrix<MT,SO2>& lhs, const DMatTSMatAddExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ResultType tmp( rhs );
      smpSchurAssign( *lhs, tmp );
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
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_FORM_VALID_MATMATADDEXPR( MT1, MT2 );
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
/*!\brief Backend implementation of the addition between a row-major dense matrix and a
//        column-major sparse matrix (\f$ A=B+C \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the addition.
// \param rhs The right-hand side sparse matrix for the addition.
// \return The sum of the two matrices.
//
// This function implements a performance optimized treatment of the addition between a row-major
// dense matrix and a column-major sparse matrix.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side sparse matrix
        , DisableIf_t< IsZero_v<MT2> &&
                       IsSame_v< ElementType_t<MT1>, ElementType_t<MT2> > >* = nullptr >
inline const DMatTSMatAddExpr<MT1,MT2>
   dmattsmatadd( const DenseMatrix<MT1,false>& lhs, const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).columns(), "Invalid number of columns" );

   return DMatTSMatAddExpr<MT1,MT2>( *lhs, *rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the addition between a row-major dense matrix and a
//        column-major zero matrix (\f$ A=B+C \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the addition.
// \param rhs The right-hand side zero matrix for the addition.
// \return Reference to the left-hand side dense matrix.
//
// This function implements a performance optimized treatment of the addition between a row-major
// dense matrix and a column-major zero matrix. It returns a reference to the left-hand side dense
// matrix.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side sparse matrix
        , EnableIf_t< IsZero_v<MT2> &&
                      IsSame_v< ElementType_t<MT1>, ElementType_t<MT2> > >* = nullptr >
inline const MT1&
   dmattsmatadd( const DenseMatrix<MT1,false>& lhs, const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( rhs );

   BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).columns(), "Invalid number of columns" );

   return (*lhs);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition operator for the addition of a row-major dense matrix and a column-major
//        sparse matrix (\f$ A=B+C \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the matrix addition.
// \param rhs The right-hand side sparse matrix to be added to the left-hand side matrix.
// \return The sum of the two matrices.
// \exception std::invalid_argument Matrix sizes do not match.
//
// This operator represents the addition of a row-major dense matrix and a column-major sparse
// matrix:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,rowMajor> A, C;
   blaze::CompressedMatrix<double,columnMajor> B;
   // ... Resizing and initialization
   C = A + B;
   \endcode

// The operator returns an expression representing a dense matrix of the higher-order element
// type of the two involved matrix element types \a MT1::ElementType and \a MT2::ElementType.
// Both matrix types \a MT1 and \a MT2 as well as the two element types \a MT1::ElementType
// and \a MT2::ElementType have to be supported by the AddTrait class template.\n
// In case the current sizes of the two given matrices don't match, a \a std::invalid_argument
// is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
inline decltype(auto)
   operator+( const DenseMatrix<MT1,false>& lhs, const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   if( (*lhs).rows() != (*rhs).rows() || (*lhs).columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   return dmattsmatadd( *lhs, *rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition operator for the addition of a column-major sparse matrix and a row-major
//        dense matrix (\f$ A=B+C \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side sparse matrix for the matrix addition.
// \param rhs The right-hand side dense matrix to be added to the left-hand side matrix.
// \return The sum of the two matrices.
// \exception std::invalid_argument Matrix sizes do not match.
//
// This operator represents the addition of a column-major sparse matrix and a row-major dense
// matrix:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   blaze::CompressedMatrix<double,columnMajor> A;
   blaze::DynamicMatrix<double,rowMajor> B, C;
   // ... Resizing and initialization
   C = A + B;
   \endcode

// The operator returns an expression representing a dense matrix of the higher-order element
// type of the two involved matrix element types \a MT1::ElementType and \a MT2::ElementType.
// Both matrix types \a MT1 and \a MT2 as well as the two element types \a MT1::ElementType
// and \a MT2::ElementType have to be supported by the AddTrait class template.\n
// In case the current sizes of the two given matrices don't match, a \a std::invalid_argument
// is thrown.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
inline decltype(auto)
   operator+( const SparseMatrix<MT1,true>& lhs, const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   if( (*lhs).rows() != (*rhs).rows() || (*lhs).columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   return dmattsmatadd( *rhs, *lhs );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition operator for the addition of a dense matrix-transpose sparse matrix
//        addition expression and a dense matrix (\f$ A=(B+C)+D \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix-transpose sparse matrix addition.
// \param rhs The right-hand side dense matrix.
// \return The sum of the two matrices.
//
// This operator implements a performance optimized treatment of the addition of a dense
// matrix-transpose sparse matrix addition expression to a dense matrix.
*/
template< typename MT1  // Type of the dense matrix of the left-hand side expression
        , typename MT2  // Type of the sparse matrix of the left-hand side expression
        , typename MT3  // Type of the right-hand side dense matrix
        , bool SO >     // Storage order of the right-hand side dense matrix
inline decltype(auto)
   operator+( const DMatTSMatAddExpr<MT1,MT2>& lhs, const DenseMatrix<MT3,SO>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return ( lhs.leftOperand() + (*rhs) ) + lhs.rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction operator for the subtraction of a dense matrix-transpose sparse matrix
//        addition expression and a dense matrix (\f$ A=(B+C)-D \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix-transpose sparse matrix addition.
// \param rhs The right-hand side dense matrix.
// \return The difference of the two matrices.
//
// This operator implements a performance optimized treatment of the subtraction of a dense
// matrix-transpose sparse matrix addition expression and a dense matrix.
*/
template< typename MT1  // Type of the dense matrix of the left-hand side expression
        , typename MT2  // Type of the sparse matrix of the left-hand side expression
        , typename MT3  // Type of the right-hand side dense matrix
        , bool SO >     // Storage order of the right-hand side dense matrix
inline decltype(auto)
   operator-( const DMatTSMatAddExpr<MT1,MT2>& lhs, const DenseMatrix<MT3,SO>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return ( lhs.leftOperand() - (*rhs) ) + lhs.rightOperand();
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
