//=================================================================================================
/*!
//  \file blaze/math/expressions/DMatTDMatSchurExpr.h
//  \brief Header file for the dense matrix/transpose dense matrix Schur product expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DMATTDMATSCHUREXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DMATTDMATSCHUREXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/Identity.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
#include <blaze/math/constraints/SchurExpr.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/constraints/Zero.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/SchurExpr.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/sparse/Forward.h>
#include <blaze/math/traits/SchurTrait.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsCommutative.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsOperation.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/typetraits/IsTemporary.h>
#include <blaze/math/typetraits/IsUniLower.h>
#include <blaze/math/typetraits/IsUniUpper.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/system/Blocking.h>
#include <blaze/system/Thresholds.h>
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
//  CLASS DMATTDMATSCHUREXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for dense matrix-transpose dense matrix Schur product.
// \ingroup dense_matrix_expression
//
// The DMatTDMatSchurExpr class represents the compile time expression for Schur product between
// a row-major dense matrix and a column-major dense matrix.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
class DMatTDMatSchurExpr
   : public SchurExpr< DenseMatrix< DMatTDMatSchurExpr<MT1,MT2>, false > >
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   using RT1 = ResultType_t<MT1>;     //!< Result type of the left-hand side dense matrix expression.
   using RT2 = ResultType_t<MT2>;     //!< Result type of the right-hand side dense matrix expression.
   using RN1 = ReturnType_t<MT1>;     //!< Return type of the left-hand side dense matrix expression.
   using RN2 = ReturnType_t<MT2>;     //!< Return type of the right-hand side dense matrix expression.
   using CT1 = CompositeType_t<MT1>;  //!< Composite type of the left-hand side dense matrix expression.
   using CT2 = CompositeType_t<MT2>;  //!< Composite type of the right-hand side dense matrix expression.
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
   //! Compilation switch for the serial evaluation strategy of the Schur product expression.
   /*! The \a useAssign compile time constant expression represents a compilation switch for the
       serial evaluation strategy of the Schur product expression. In case either of the two dense
       matrix operands requires an intermediate evaluation or the subscript operator can only
       return by value, \a useAssign will be set to 1 and the Schur product expression will be
       evaluated via the \a assign function family. Otherwise \a useAssign will be set to 0
       and the expression will be evaluated via the function call operator. */
   static constexpr bool useAssign =
      ( RequiresEvaluation_v<MT1> || RequiresEvaluation_v<MT2> || !returnExpr );

   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT >
   static constexpr bool UseAssign_v = useAssign;
   /*! \endcond */
   //**********************************************************************************************

   //**Parallel evaluation strategy****************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! This variable template is a helper for the selection of the parallel evaluation strategy.
       In case at least one of the two matrix operands is not SMP assignable and at least one of
       the two operands requires an intermediate evaluation, the variable is set to 1 and the
       expression specific evaluation strategy is selected. Otherwise the variable is set to 0
       and the default strategy is chosen. */
   template< typename MT >
   static constexpr bool UseSMPAssign_v =
      ( ( !MT1::smpAssignable || !MT2::smpAssignable ) && useAssign );
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this DMatTDMatSchurExpr instance.
   using This = DMatTDMatSchurExpr<MT1,MT2>;

   //! Base type of this DMatTDMatSchurExpr instance.
   using BaseType = SchurExpr< DenseMatrix<This,false> >;

   using ResultType    = SchurTrait_t<RT1,RT2>;        //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<ResultType>;   //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<ResultType>;    //!< Resulting element type.

   //! Return type for expression template evaluations.
   using ReturnType = const If_t< returnExpr, ExprReturnType, ElementType >;

   //! Data type for composite expression templates.
   using CompositeType = const ResultType;

   //! Composite type of the left-hand side dense matrix expression.
   using LeftOperand = If_t< IsExpression_v<MT1>, const MT1, const MT1& >;

   //! Composite type of the right-hand side dense matrix expression.
   using RightOperand = If_t< IsExpression_v<MT2>, const MT2, const MT2& >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = false;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = ( MT1::smpAssignable && MT2::smpAssignable );
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DMatTDMatSchurExpr class.
   //
   // \param lhs The left-hand side operand of the Schur product expression.
   // \param rhs The right-hand side operand of the Schur product expression.
   */
   inline DMatTDMatSchurExpr( const MT1& lhs, const MT2& rhs ) noexcept
      : lhs_( lhs )  // Left-hand side dense matrix of the Schur product expression
      , rhs_( rhs )  // Right-hand side dense matrix of the Schur product expression
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

   //**Left operand access*************************************************************************
   /*!\brief Returns the left-hand side dense matrix operand.
   //
   // \return The left-hand side dense matrix operand.
   */
   inline LeftOperand leftOperand() const noexcept {
      return lhs_;
   }
   //**********************************************************************************************

   //**Right operand access************************************************************************
   /*!\brief Returns the right-hand side transpose dense matrix operand.
   //
   // \return The right-hand side transpose dense matrix operand.
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
      return ( IsExpression_v<MT1> && ( RequiresEvaluation_v<MT1> ? lhs_.isAliased( alias ) : lhs_.canAlias( alias ) ) ) ||
             ( IsExpression_v<MT2> && ( RequiresEvaluation_v<MT2> ? rhs_.isAliased( alias ) : rhs_.canAlias( alias ) ) );
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

   //**********************************************************************************************
   /*!\brief Returns whether the operands of the expression are properly aligned in memory.
   //
   // \return \a true in case the operands are aligned, \a false if not.
   */
   inline bool isAligned() const noexcept {
      return lhs_.isAligned() && rhs_.isAligned();
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const noexcept {
      return lhs_.canSMPAssign() || rhs_.canSMPAssign() ||
             ( rows() * columns() >= SMP_DMATTDMATSCHUR_THRESHOLD );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  lhs_;  //!< Left-hand side dense matrix of the Schur product expression.
   RightOperand rhs_;  //!< Right-hand side dense matrix of the Schur product expression.
   //**********************************************************************************************

   //**Assignment to dense matrices****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense matrix-transpose dense matrix Schur product to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Schur product expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix-transpose
   // dense matrix Schur product expression to a dense matrix. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case neither
   // of the two operands requires an intermediate evaluation.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto assign( DenseMatrix<MT,SO>& lhs, const DMatTDMatSchurExpr& rhs )
      -> DisableIf_t< UseAssign_v<MT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      constexpr size_t block( BLOCK_SIZE );

      const size_t m( rhs.rows() );
      const size_t n( rhs.columns() );

      for( size_t ii=0UL; ii<m; ii+=block ) {
         const size_t iend( ( m < ii+block )?( m ):( ii+block ) );
         for( size_t jj=0UL; jj<n; jj+=block ) {
            const size_t jend( ( n < jj+block )?( n ):( jj+block ) );
            for( size_t i=ii; i<iend; ++i ) {
               for( size_t j=jj; j<jend; ++j ) {
                  (*lhs)(i,j) = rhs.lhs_(i,j) * rhs.rhs_(i,j);
               }
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to dense matrices****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a non-commutative dense matrix-transpose dense matrix Schur product
   //        to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Schur product expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a non-commutative dense
   // matrix-transpose dense matrix Schur product expression to a dense matrix. Due to the
   // explicit application of the SFINAE principle, this function can only be selected by
   // the compiler in case either of the two operands requires an intermediate evaluation.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto assign( DenseMatrix<MT,SO>& lhs, const DMatTDMatSchurExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT> && !IsCommutative_v<MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( !IsOperation_v<MT1> && isSame( *lhs, rhs.lhs_ ) ) {
         schurAssign( *lhs, rhs.rhs_ );
      }
      else {
         CT1 A( serial( rhs.lhs_ ) );
         CT2 B( serial( rhs.rhs_ ) );
         assign( *lhs, A % B );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to dense matrices****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a commutative dense matrix-transpose dense matrix Schur product to
   //        a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Schur product expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a commutative dense
   // matrix-transpose dense matrix Schur product expression to a dense matrix. Due to the
   // explicit application of the SFINAE principle, this function can only be selected by
   // the compiler in case either of the two operands requires an intermediate evaluation.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto assign( DenseMatrix<MT,SO>& lhs, const DMatTDMatSchurExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT> && IsCommutative_v<MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( !IsOperation_v<MT1> && isSame( *lhs, rhs.lhs_ ) ) {
         schurAssign( *lhs, rhs.rhs_ );
      }
      else if( !IsOperation_v<MT2> && isSame( *lhs, rhs.rhs_ ) ) {
         schurAssign( *lhs, rhs.lhs_ );
      }
      else if( !RequiresEvaluation_v<MT2> ) {
         assign     ( *lhs, rhs.rhs_ );
         schurAssign( *lhs, rhs.lhs_ );
      }
      else {
         assign     ( *lhs, rhs.lhs_ );
         schurAssign( *lhs, rhs.rhs_ );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse matrices***************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense matrix-transpose dense matrix Schur product to a sparse matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side Schur product expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix-transpose
   // dense matrix Schur product expression to a sparse matrix.
   */
   template< typename MT  // Type of the target sparse matrix
           , bool SO >    // Storage order of the target sparse matrix
   friend inline void assign( SparseMatrix<MT,SO>& lhs, const DMatTDMatSchurExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      using TmpType = If_t< SO, OppositeType, ResultType >;

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
   /*!\brief Addition assignment of a dense matrix-transpose dense matrix Schur product to a
   //        dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Schur product expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense matrix-
   // transpose dense matrix Schur product expression to a dense matrix. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler
   // in case neither of the two operands requires an intermediate evaluation.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto addAssign( DenseMatrix<MT,SO>& lhs, const DMatTDMatSchurExpr& rhs )
      -> DisableIf_t< UseAssign_v<MT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      constexpr size_t block( BLOCK_SIZE );

      const size_t m( rhs.rows() );
      const size_t n( rhs.columns() );

      for( size_t ii=0UL; ii<m; ii+=block ) {
         const size_t iend( ( m < ii+block )?( m ):( ii+block ) );
         for( size_t jj=0UL; jj<n; jj+=block ) {
            const size_t jend( ( n < jj+block )?( n ):( jj+block ) );
            for( size_t i=ii; i<iend; ++i ) {
               for( size_t j=jj; j<jend; ++j ) {
                  (*lhs)(i,j) += rhs.lhs_(i,j) * rhs.rhs_(i,j);
               }
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense matrices*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a dense matrix-transpose dense matrix Schur product to a
   //        dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Schur product expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense matrix-
   // transpose dense matrix Schur product expression to a dense matrix. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler
   // in case either of the two operands requires an intermediate evaluation.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto addAssign( DenseMatrix<MT,SO>& lhs, const DMatTDMatSchurExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ResultType tmp( serial( rhs ) );
      addAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse matrices******************************************************
   // No special implementation for the addition assignment to sparse matrices.
   //**********************************************************************************************

   //**Subtraction assignment to dense matrices****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a dense matrix-transpose dense matrix Schur product to a
   //        dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Schur product expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a dense
   // matrix-transpose dense matrix Schur product expression to a dense matrix. Due to the
   // explicit application of the SFINAE principle, this function can only be selected by
   // the compiler in case neither of the two operands requires an intermediate evaluation.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto subAssign( DenseMatrix<MT,SO>& lhs, const DMatTDMatSchurExpr& rhs )
      -> DisableIf_t< UseAssign_v<MT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      constexpr size_t block( BLOCK_SIZE );

      const size_t m( rhs.rows() );
      const size_t n( rhs.columns() );

      for( size_t ii=0UL; ii<m; ii+=block ) {
         const size_t iend( ( m < ii+block )?( m ):( ii+block ) );
         for( size_t jj=0UL; jj<n; jj+=block ) {
            const size_t jend( ( n < jj+block )?( n ):( jj+block ) );
            for( size_t i=ii; i<iend; ++i ) {
               for( size_t j=jj; j<jend; ++j ) {
                  (*lhs)(i,j) -= rhs.lhs_(i,j) * rhs.rhs_(i,j);
               }
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to dense matrices****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a dense matrix-transpose dense matrix Schur product to a
   //        dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Schur product expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a dense
   // matrix-transpose dense matrix Schur product expression to a dense matrix. Due to the
   // explicit application of the SFINAE principle, this function can only be selected by
   // the compiler in case either of the two operands requires an intermediate evaluation.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto subAssign( DenseMatrix<MT,SO>& lhs, const DMatTDMatSchurExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ResultType tmp( serial( rhs ) );
      subAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse matrices***************************************************
   // No special implementation for the subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**Schur product assignment to dense matrices**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Schur product assignment of a dense matrix-transpose dense matrix Schur product to
   //        a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Schur product expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a dense
   // matrix-transpose dense matrix Schur product expression to a dense matrix. Due to the
   // explicit application of the SFINAE principle, this function can only be selected by
   // the compiler in case neither of the two operands requires an intermediate evaluation.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto schurAssign( DenseMatrix<MT,SO>& lhs, const DMatTDMatSchurExpr& rhs )
      -> DisableIf_t< UseAssign_v<MT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      constexpr size_t block( BLOCK_SIZE );

      const size_t m( rhs.rows() );
      const size_t n( rhs.columns() );

      for( size_t ii=0UL; ii<m; ii+=block ) {
         const size_t iend( ( m < ii+block )?( m ):( ii+block ) );
         for( size_t jj=0UL; jj<n; jj+=block ) {
            const size_t jend( ( n < jj+block )?( n ):( jj+block ) );
            for( size_t i=ii; i<iend; ++i ) {
               for( size_t j=jj; j<jend; ++j ) {
                  (*lhs)(i,j) *= rhs.lhs_(i,j) * rhs.rhs_(i,j);
               }
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Schur product assignment to dense matrices**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Schur product assignment of a non-commutative dense matrix-transpose dense matrix
   //        Schur product to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Schur product expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a
   // non-commutative dense matrix-transpose dense matrix Schur product expression to a dense
   // matrix. Due to the explicit application of the SFINAE principle, this function can only
   // be selected by the compiler in case either of the two operands requires an intermediate
   // evaluation.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto schurAssign( DenseMatrix<MT,SO>& lhs, const DMatTDMatSchurExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT> && !IsCommutative_v<MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( ResultType, SO );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ResultType tmp( serial( rhs ) );
      schurAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Schur product assignment to dense matrices**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Schur product assignment of a commutative dense matrix-transpose dense matrix
   //        Schur product to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Schur product expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a
   // commutative dense matrix-transpose dense matrix Schur product expression to a dense
   // matrix. Due to the explicit application of the SFINAE principle, this function can only
   // be selected by the compiler in case either of the two operands requires an intermediate
   // evaluation.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto schurAssign( DenseMatrix<MT,SO>& lhs, const DMatTDMatSchurExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT> && IsCommutative_v<MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( !RequiresEvaluation_v<MT2> ) {
         schurAssign( *lhs, rhs.rhs_ );
         schurAssign( *lhs, rhs.lhs_ );
      }
      else {
         schurAssign( *lhs, rhs.lhs_ );
         schurAssign( *lhs, rhs.rhs_ );
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

   //**SMP assignment to dense matrices************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a non-commutative dense matrix-transpose dense matrix Schur product
   //        to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Schur product expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a non-commutative
   // dense matrix-transpose dense matrix Schur product expression to a dense matrix. Due to
   // the explicit application of the SFINAE principle, this function can only be selected by
   // the compiler in case the expression specific parallel evaluation strategy is selected.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto smpAssign( DenseMatrix<MT,SO>& lhs, const DMatTDMatSchurExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT> && !IsCommutative_v<MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( !IsOperation_v<MT1> && isSame( *lhs, rhs.lhs_ ) ) {
         smpSchurAssign( *lhs, rhs.rhs_ );
      }
      else {
         CT1 A( rhs.lhs_ );
         CT2 B( rhs.rhs_ );
         smpAssign( *lhs, A % B );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to dense matrices************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a commutative dense matrix-transpose dense matrix Schur product
   //        to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Schur product expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a commutative
   // dense matrix-transpose dense matrix Schur product expression to a dense matrix. Due to
   // the explicit application of the SFINAE principle, this function can only be selected by
   // the compiler in case the expression specific parallel evaluation strategy is selected.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto smpAssign( DenseMatrix<MT,SO>& lhs, const DMatTDMatSchurExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT> && IsCommutative_v<MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( !IsOperation_v<MT1> && isSame( *lhs, rhs.lhs_ ) ) {
         smpSchurAssign( *lhs, rhs.rhs_ );
      }
      else if( !IsOperation_v<MT2> && isSame( *lhs, rhs.rhs_ ) ) {
         smpSchurAssign( *lhs, rhs.lhs_ );
      }
      else if( !RequiresEvaluation_v<MT2> ) {
         smpAssign     ( *lhs, rhs.rhs_ );
         smpSchurAssign( *lhs, rhs.lhs_ );
      }
      else {
         smpAssign     ( *lhs, rhs.lhs_ );
         smpSchurAssign( *lhs, rhs.rhs_ );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to sparse matrices***********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a dense matrix-transpose dense matrix Schur product to a sparse
   //        matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side Schur product expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a dense matrix-
   // transpose dense matrix Schur product expression to a sparse matrix. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler
   // in case the expression specific parallel evaluation strategy is selected.
   */
   template< typename MT  // Type of the target sparse matrix
           , bool SO >    // Storage order of the target sparse matrix
   friend inline auto smpAssign( SparseMatrix<MT,SO>& lhs, const DMatTDMatSchurExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT> >
   {
      BLAZE_FUNCTION_TRACE;

      using TmpType = If_t< SO, OppositeType, ResultType >;

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
   /*!\brief SMP addition assignment of a dense matrix-transpose dense matrix Schur product to
   //        a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Schur product expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a dense
   // matrix-transpose dense matrix Schur product expression to a dense matrix. Due to the
   // explicit application of the SFINAE principle, this function can only be selected by
   // the compiler in case the expression specific parallel evaluation strategy is selected.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto smpAddAssign( DenseMatrix<MT,SO>& lhs, const DMatTDMatSchurExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ResultType tmp( rhs );
      smpAddAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to sparse matrices**************************************************
   // No special implementation for the SMP addition assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP subtraction assignment to dense matrices************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a dense matrix-transpose dense matrix Schur product
   //        to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Schur product expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a dense
   // matrix-transpose dense matrix Schur product expression to a dense matrix. Due to the
   // explicit application of the SFINAE principle, this function can only be selected by the
   // compiler in case the expression specific parallel evaluation strategy is selected.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto smpSubAssign( DenseMatrix<MT,SO>& lhs, const DMatTDMatSchurExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ResultType tmp( rhs );
      smpSubAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to sparse matrices***********************************************
   // No special implementation for the SMP subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP Schur product assignment to dense matrices**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP Schur product assignment of a non-commutative dense matrix-transpose dense matrix
   //        Schur product to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Schur product expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized SMP Schur product assignment of a
   // non-commutative dense matrix-transpose dense matrix Schur product expression to a dense
   // matrix. Due to the explicit application of the SFINAE principle, this function can only
   // be selected by the compiler in case the expression specific parallel evaluation strategy
   // is selected.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto smpSchurAssign( DenseMatrix<MT,SO>& lhs, const DMatTDMatSchurExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT> && !IsCommutative_v<MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( ResultType, SO );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ResultType tmp( rhs );
      smpSchurAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP Schur product assignment to dense matrices**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP Schur product assignment of a commutative dense matrix-transpose dense matrix
   //        Schur product to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Schur product expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized SMP Schur product assignment of a
   // commutative dense matrix-transpose dense matrix Schur product expression to a dense
   // matrix. Due to the explicit application of the SFINAE principle, this function can only
   // be selected by the compiler in case the expression specific parallel evaluation strategy
   // is selected.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto smpSchurAssign( DenseMatrix<MT,SO>& lhs, const DMatTDMatSchurExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT> && IsCommutative_v<MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( !RequiresEvaluation_v<MT2> ) {
         smpSchurAssign( *lhs, rhs.rhs_ );
         smpSchurAssign( *lhs, rhs.lhs_ );
      }
      else {
         smpSchurAssign( *lhs, rhs.lhs_ );
         smpSchurAssign( *lhs, rhs.rhs_ );
      }
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
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT2 );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER( MT1, MT2 );
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
/*!\brief Backend implementation of the Schur product between a row-major dense matrix and
//        a column-major dense matrix (\f$ A=B \circ C \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the Schur product.
// \param rhs The right-hand side sparse matrix for the Schur product.
// \return The Schur product of the two matrices.
//
// This function implements a performance optimized treatment of the Schur product between a
// row-major dense matrix and a column-major dense matrix.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , EnableIf_t< !IsSymmetric_v<MT1> &&
                      !IsSymmetric_v<MT2> &&
                      !( IsUniLower_v<MT1> && IsUniUpper_v<MT2> ) &&
                      !( IsUniUpper_v<MT1> && IsUniLower_v<MT2> ) &&
                      !( IsStrictlyLower_v<MT1> && IsUpper_v<MT2> ) &&
                      !( IsStrictlyUpper_v<MT1> && IsLower_v<MT2> ) &&
                      !( IsLower_v<MT1> && IsStrictlyUpper_v<MT2> ) &&
                      !( IsUpper_v<MT1> && IsStrictlyLower_v<MT2> ) >* = nullptr >
inline const DMatTDMatSchurExpr<MT1,MT2>
   dmattdmatschur( const DenseMatrix<MT1,false>& lhs, const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).columns(), "Invalid number of columns" );

   return DMatTDMatSchurExpr<MT1,MT2>( *lhs, *rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the Schur product between a symmetric row-major dense matrix
//        and a column-major dense matrix (\f$ A=B \circ C \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the Schur product.
// \param rhs The right-hand side dense matrix for the Schur product.
// \return The Schur product of the two matrices.
//
// This function implements a performance optimized treatment of the Schur product of a symmetric
// row-major dense matrix and a column-major dense matrix.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , EnableIf_t< IsSymmetric_v<MT1> &&
                      !IsSymmetric_v<MT2> &&
                      !( IsUniLower_v<MT1> && IsUniUpper_v<MT2> ) &&
                      !( IsUniUpper_v<MT1> && IsUniLower_v<MT2> ) &&
                      !( IsStrictlyLower_v<MT1> && IsUpper_v<MT2> ) &&
                      !( IsStrictlyUpper_v<MT1> && IsLower_v<MT2> ) &&
                      !( IsLower_v<MT1> && IsStrictlyUpper_v<MT2> ) &&
                      !( IsUpper_v<MT1> && IsStrictlyLower_v<MT2> ) >* = nullptr >
inline decltype(auto)
   dmattdmatschur( const DenseMatrix<MT1,false>& lhs, const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).columns(), "Invalid number of columns" );

   return trans( *lhs ) % *rhs;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation for the Schur product between a row-major dense matrix and a
//        symmetric column-major dense matrix (\f$ A=B \circ C \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the Schur product.
// \param rhs The right-hand side dense matrix for the Schur product.
// \return The sum of the two matrices.
//
// This function implements a performance optimized treatment of the Schur product of a (potentially
// symmetric) row-major dense matrix and a symmetric column-major dense matrix.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , EnableIf_t< IsSymmetric_v<MT2> &&
                      !( IsUniLower_v<MT1> && IsUniUpper_v<MT2> ) &&
                      !( IsUniUpper_v<MT1> && IsUniLower_v<MT2> ) &&
                      !( IsStrictlyLower_v<MT1> && IsUpper_v<MT2> ) &&
                      !( IsStrictlyUpper_v<MT1> && IsLower_v<MT2> ) &&
                      !( IsLower_v<MT1> && IsStrictlyUpper_v<MT2> ) &&
                      !( IsUpper_v<MT1> && IsStrictlyLower_v<MT2> ) >* = nullptr >
inline decltype(auto)
   dmattdmatschur( const DenseMatrix<MT1,false>& lhs, const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).columns(), "Invalid number of columns" );

   return (*lhs) % trans( *rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the Schur product between a unitriangular row-major dense
//        matrix and a unitriangular column-major dense matrix (\f$ A=B \circ C \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the Schur product.
// \param rhs The right-hand side dense matrix for the Schur product.
// \return The resulting identity matrix.
//
// This function implements a performance optimized treatment of the Schur product between a
// unitriangular row-major dense matrix and a unitriangular column-major dense matrix. It
// returns an identity matrix.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , EnableIf_t< ( IsUniLower_v<MT1> && IsUniUpper_v<MT2> ) ||
                      ( IsUniUpper_v<MT1> && IsUniLower_v<MT2> ) >* = nullptr >
inline decltype(auto)
   dmattdmatschur( const DenseMatrix<MT1,false>& lhs, const DenseMatrix<MT2,true>& rhs )
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
//        dense matrix and a (strictly) triangular column-major dense matrix (\f$ A=B \circ C \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the Schur product.
// \param rhs The right-hand side dense matrix for the Schur product.
// \return The resulting zero matrix.
//
// This function implements a performance optimized treatment of the Schur product between a
// (strictly) triangular row-major dense matrix and a (strictly) triangular column-major dense
// matrix. It returns a zero matrix.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , EnableIf_t< ( IsStrictlyLower_v<MT1> && IsUpper_v<MT2> ) ||
                      ( IsStrictlyUpper_v<MT1> && IsLower_v<MT2> ) ||
                      ( IsLower_v<MT1> && IsStrictlyUpper_v<MT2> ) ||
                      ( IsUpper_v<MT1> && IsStrictlyLower_v<MT2> ) >* = nullptr >
inline decltype(auto)
   dmattdmatschur( const DenseMatrix<MT1,false>& lhs, const DenseMatrix<MT2,true>& rhs )
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
/*!\brief Operator for the Schur product of a row-major and a column-major dense matrix
//        (\f$ A=B \circ C \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the Schur product.
// \param rhs The right-hand side dense matrix for the Schur product.
// \return The Schur product of the two matrices.
// \exception std::invalid_argument Matrix sizes do not match.
//
// This operator represents the Schur product of a row-major and a column-major dense matrix:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,rowMajor> A, C;
   blaze::DynamicMatrix<double,columnMajor> B;
   // ... Resizing and initialization
   C = A % B;
   \endcode

// The operator returns an expression representing a dense matrix of the higher-order element
// type of the two involved matrix element types \a MT1::ElementType and \a MT2::ElementType.
// Both matrix types \a MT1 and \a MT2 as well as the two element types \a MT1::ElementType
// and \a MT2::ElementType have to be supported by the MultTrait class template.\n
// In case the current number of rows and columns of the two given  matrices don't match, a
// \a std::invalid_argument is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
inline decltype(auto)
   operator%( const DenseMatrix<MT1,false>& lhs, const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   if( (*lhs).rows() != (*rhs).rows() || (*lhs).columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   return dmattdmatschur( *lhs, *rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the Schur product between a column-major dense matrix and
//        a row-major dense matrix (\f$ A=B \circ C \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the Schur product.
// \param rhs The right-hand side sparse matrix for the Schur product.
// \return The Schur product of the two matrices.
//
// This function implements a performance optimized treatment of the Schur product between a
// column-major dense matrix and a row-major dense matrix.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , EnableIf_t< !IsSymmetric_v<MT1> &&
                      !IsSymmetric_v<MT2> &&
                      !( IsUniLower_v<MT1> && IsUniUpper_v<MT2> ) &&
                      !( IsUniUpper_v<MT1> && IsUniLower_v<MT2> ) &&
                      !( IsStrictlyLower_v<MT1> && IsUpper_v<MT2> ) &&
                      !( IsStrictlyUpper_v<MT1> && IsLower_v<MT2> ) &&
                      !( IsLower_v<MT1> && IsStrictlyUpper_v<MT2> ) &&
                      !( IsUpper_v<MT1> && IsStrictlyLower_v<MT2> ) >* = nullptr >
inline const DMatTDMatSchurExpr<MT1,MT2>
   tdmatdmatschur( const DenseMatrix<MT1,true>& lhs, const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).columns(), "Invalid number of columns" );

   return DMatTDMatSchurExpr<MT1,MT2>( *lhs, *rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the Schur product between a column-major dense matrix and
//        a symmetric row-major dense matrix (\f$ A=B \circ C \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the Schur product.
// \param rhs The right-hand side dense matrix for the Schur product.
// \return The Schur product of the two matrices.
//
// This function implements a performance optimized treatment of the Schur product of a
// column-major dense matrix and a symmetric row-major dense matrix.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , EnableIf_t< !IsSymmetric_v<MT1> &&
                      IsSymmetric_v<MT2> &&
                      !( IsUniLower_v<MT1> && IsUniUpper_v<MT2> ) &&
                      !( IsUniUpper_v<MT1> && IsUniLower_v<MT2> ) &&
                      !( IsStrictlyLower_v<MT1> && IsUpper_v<MT2> ) &&
                      !( IsStrictlyUpper_v<MT1> && IsLower_v<MT2> ) &&
                      !( IsLower_v<MT1> && IsStrictlyUpper_v<MT2> ) &&
                      !( IsUpper_v<MT1> && IsStrictlyLower_v<MT2> ) >* = nullptr >
inline decltype(auto)
   tdmatdmatschur( const DenseMatrix<MT1,true>& lhs, const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).columns(), "Invalid number of columns" );

   return (*lhs) % trans( *rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation for the Schur product between a symmetric column-major dense
//        matrix and a row-major dense matrix (\f$ A=B \circ C \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the Schur product.
// \param rhs The right-hand side dense matrix for the Schur product.
// \return The sum of the two matrices.
//
// This function implements a performance optimized treatment of the Schur product of a symmetric
// column-major dense matrix and a (potentially symmetric) row-major dense matrix.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , EnableIf_t< IsSymmetric_v<MT1> &&
                      !( IsUniLower_v<MT1> && IsUniUpper_v<MT2> ) &&
                      !( IsUniUpper_v<MT1> && IsUniLower_v<MT2> ) &&
                      !( IsStrictlyLower_v<MT1> && IsUpper_v<MT2> ) &&
                      !( IsStrictlyUpper_v<MT1> && IsLower_v<MT2> ) &&
                      !( IsLower_v<MT1> && IsStrictlyUpper_v<MT2> ) &&
                      !( IsUpper_v<MT1> && IsStrictlyLower_v<MT2> ) >* = nullptr >
inline decltype(auto)
   tdmatdmatschur( const DenseMatrix<MT1,true>& lhs, const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).columns(), "Invalid number of columns" );

   return trans( *lhs ) % (*rhs);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the Schur product between a unitriangular column-major dense
//        matrix and a unitriangular row-major dense matrix (\f$ A=B \circ C \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the Schur product.
// \param rhs The right-hand side dense matrix for the Schur product.
// \return The resulting identity matrix.
//
// This function implements a performance optimized treatment of the Schur product between a
// unitriangular column-major dense matrix and a unitriangular row-major dense matrix. It
// returns an identity matrix.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , EnableIf_t< ( IsUniLower_v<MT1> && IsUniUpper_v<MT2> ) ||
                      ( IsUniUpper_v<MT1> && IsUniLower_v<MT2> ) >* = nullptr >
inline decltype(auto)
   tdmatdmatschur( const DenseMatrix<MT1,true>& lhs, const DenseMatrix<MT2,false>& rhs )
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
/*!\brief Backend implementation of the Schur product between a (strictly) triangular column-major
//        dense matrix and a (strictly) triangular row-major dense matrix (\f$ A=B \circ C \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the Schur product.
// \param rhs The right-hand side dense matrix for the Schur product.
// \return The resulting zero matrix.
//
// This function implements a performance optimized treatment of the Schur product between a
// (strictly) triangular column-major dense matrix and a (strictly) triangular row-major dense
// matrix. It returns a zero matrix.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , EnableIf_t< ( IsStrictlyLower_v<MT1> && IsUpper_v<MT2> ) ||
                      ( IsStrictlyUpper_v<MT1> && IsLower_v<MT2> ) ||
                      ( IsLower_v<MT1> && IsStrictlyUpper_v<MT2> ) ||
                      ( IsUpper_v<MT1> && IsStrictlyLower_v<MT2> ) >* = nullptr >
inline decltype(auto)
   tdmatdmatschur( const DenseMatrix<MT1,true>& lhs, const DenseMatrix<MT2,false>& rhs )
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
/*!\brief Operator for the Schur product of a column-major and a row-major dense matrix
//        (\f$ A=B \circ C \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the Schur product.
// \param rhs The right-hand side dense matrix for the Schur product.
// \return The Schur product of the two matrices.
// \exception std::invalid_argument Matrix sizes do not match
//
// This operator represents the Schur product of a column-major and a row-major dense matrix:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> A;
   blaze::DynamicMatrix<double,rowMajor> B, C;
   // ... Resizing and initialization
   C = A % B;
   \endcode

// The operator returns an expression representing a dense matrix of the higher-order element
// type of the two involved matrix element types \a MT1::ElementType and \a MT2::ElementType.
// Both matrix types \a MT1 and \a MT2 as well as the two element types \a MT1::ElementType
// and \a MT2::ElementType have to be supported by the MultTrait class template.\n
// In case the current number of rows and columns of the two given  matrices don't match, a
// \a std::invalid_argument is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
inline decltype(auto)
   operator%( const DenseMatrix<MT1,true>& lhs, const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   if( (*lhs).rows() != (*rhs).rows() || (*lhs).columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   return tdmatdmatschur( *lhs, *rhs );
}
//*************************************************************************************************




//=================================================================================================
//
//  ISALIGNED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2 >
struct IsAligned< DMatTDMatSchurExpr<MT1,MT2> >
   : public BoolConstant< IsAligned_v<MT1> && IsAligned_v<MT2> >
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
