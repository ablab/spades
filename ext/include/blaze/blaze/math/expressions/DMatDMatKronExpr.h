//=================================================================================================
/*!
//  \file blaze/math/expressions/DMatDMatKronExpr.h
//  \brief Header file for the dense matrix/dense matrix Kronecker product expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DMATDMATKRONEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DMATDMATKRONEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/math/Aliases.h>
#include <blaze/math/AlignmentFlag.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/MatMatKronExpr.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/MatMatKronExpr.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/traits/KronTrait.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsTemporary.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/util/Assert.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DMATDMATKRONEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for dense matrix-dense matrix Kronecker products.
// \ingroup dense_matrix_expression
//
// The DMatDMatKronExpr class represents the compile time expression for Kronecker products
// between dense matrices of any storage order.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO >     // Storage order
class DMatDMatKronExpr
   : public MatMatKronExpr< DenseMatrix< DMatDMatKronExpr<MT1,MT2,SO>, SO > >
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
   using ET1 = ElementType_t<MT1>;    //!< Element type of the left-hand side dense matrix expression.
   using ET2 = ElementType_t<MT2>;    //!< Element type of the right-hand side dense matrix expression.
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
   //! Type of this DMatDMatKronExpr instance.
   using This = DMatDMatKronExpr<MT1,MT2,SO>;

   //! Base type of this DMatDMatKronExpr instance.
   using BaseType = MatMatKronExpr< DenseMatrix<This,SO> >;

   using ResultType    = KronTrait_t<RT1,RT2>;         //!< Result type for expression template evaluations.
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
   static constexpr bool smpAssignable = false;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DMatDMatKronExpr class.
   //
   // \param lhs The left-hand side operand of the Kronecker product expression.
   // \param rhs The right-hand side operand of the Kronecker product expression.
   */
   inline DMatDMatKronExpr( const MT1& lhs, const MT2& rhs ) noexcept
      : lhs_( lhs )  // Left-hand side dense matrix of the Kronecker product expression
      , rhs_( rhs )  // Right-hand side dense matrix of the Kronecker product expression
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
   /*!\brief Returns the right-hand side dense matrix operand.
   //
   // \return The right-hand side dense matrix operand.
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

   //**********************************************************************************************
   /*!\brief Returns whether the operands of the expression are properly aligned in memory.
   //
   // \return \a true in case the operands are aligned, \a false if not.
   */
   inline bool isAligned() const noexcept {
      return rhs_.isAligned();
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const noexcept {
      return false;
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  lhs_;  //!< Left-hand side dense matrix of the Kronecker product expression.
   RightOperand rhs_;  //!< Right-hand side dense matrix of the Kronecker product expression.
   //**********************************************************************************************

   //**Assignment to row-major dense matrices******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense matrix-dense matrix Kronecker product to a row-major dense
   //        matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Kronecker product expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix-dense
   // matrix Kronecker product expression to a row-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline void assign( DenseMatrix<MT,false>& lhs, const DMatDMatKronExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( rhs.rows() == 0UL || rhs.columns() == 0UL ) {
         return;
      }

      CT1 A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense matrix operand
      CT2 B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side dense matrix operand

      const size_t M( B.rows()    );
      const size_t N( B.columns() );

      if( SO )
      {
         for( size_t i=0UL; i<A.rows(); ++i ) {
            for( size_t j=0UL; j<A.columns(); ++j ) {
               if( !isDefault<strict>( A(i,j) ) ) {
                  for( size_t l=0UL; l<N; ++l )
                     for( size_t k=0UL; k<M; ++k )
                        (*lhs)(i*M+k,j*N+l) = A(i,j) * B(k,l);
               }
               else {
                  for( size_t l=0UL; l<N; ++l )
                     for( size_t k=0UL; k<M; ++k )
                        reset( (*lhs)(i*M+k,j*N+l) );
               }
            }
         }
      }
      else
      {
         for( size_t i=0UL; i<A.rows(); ++i ) {
            for( size_t k=0UL; k<M; ++k ) {
               for( size_t j=0UL; j<A.columns(); ++j ) {
                  if( !isDefault<strict>( A(i,j) ) ) {
                     for( size_t l=0UL; l<N; ++l )
                        (*lhs)(i*M+k,j*N+l) = A(i,j) * B(k,l);
                  }
                  else {
                     for( size_t l=0UL; l<N; ++l )
                        reset( (*lhs)(i*M+k,j*N+l) );
                  }
               }
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to column-major dense matrices***************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense matrix-dense matrix Kronecker product to a column-major dense
   //        matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Kronecker product expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix-dense
   // matrix Kronecker product expression to a column-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline void assign( DenseMatrix<MT,true>& lhs, const DMatDMatKronExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( rhs.rows() == 0UL || rhs.columns() == 0UL ) {
         return;
      }

      CT1 A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense matrix operand
      CT2 B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side dense matrix operand

      const size_t M( B.rows()    );
      const size_t N( B.columns() );

      if( SO )
      {
         for( size_t j=0UL; j<A.columns(); ++j ) {
            for( size_t l=0UL; l<N; ++l ) {
               for( size_t i=0UL; i<A.rows(); ++i ) {
                  if( !isDefault<strict>( A(i,j) ) ) {
                     for( size_t k=0UL; k<M; ++k )
                        (*lhs)(i*M+k,j*N+l) = A(i,j) * B(k,l);
                  }
                  else {
                     for( size_t k=0UL; k<M; ++k )
                        reset( (*lhs)(i*M+k,j*N+l) );
                  }
               }
            }
         }
      }
      else
      {
         for( size_t j=0UL; j<A.columns(); ++j ) {
            for( size_t i=0UL; i<A.rows(); ++i ) {
               if( !isDefault<strict>( A(i,j) ) ) {
                  for( size_t k=0UL; k<M; ++k )
                     for( size_t l=0UL; l<N; ++l )
                        (*lhs)(i*M+k,j*N+l) = A(i,j) * B(k,l);
               }
               else {
                  for( size_t k=0UL; k<M; ++k )
                     for( size_t l=0UL; l<N; ++l )
                        reset( (*lhs)(i*M+k,j*N+l) );
               }
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse matrices***************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense matrix-dense matrix Kronecker product to a sparse matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side Kronecker product expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix-dense
   // matrix Kronecker product expression to a sparse matrix.
   */
   template< typename MT  // Type of the target sparse matrix
           , bool SO2 >   // Storage order of the target sparse matrix
   friend inline void assign( SparseMatrix<MT,SO2>& lhs, const DMatDMatKronExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      using TmpType = If_t< SO == SO2, ResultType, OppositeType >;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( OppositeType );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( ResultType, SO );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( OppositeType, !SO );
      BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER( MT, TmpType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( TmpType );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const TmpType tmp( serial( rhs ) );
      assign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to row-major dense matrices*********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a dense matrix-dense matrix Kronecker product to a row-major
   //        dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Kronecker product expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense matrix-
   // dense matrix Kronecker product expression to a row-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline void addAssign( DenseMatrix<MT,false>& lhs, const DMatDMatKronExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( rhs.rows() == 0UL || rhs.columns() == 0UL ) {
         return;
      }

      CT1 A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense matrix operand
      CT2 B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side dense matrix operand

      const size_t M( B.rows()    );
      const size_t N( B.columns() );

      if( SO )
      {
         for( size_t i=0UL; i<A.rows(); ++i )
            for( size_t j=0UL; j<A.columns(); ++j )
               if( !isDefault<strict>( A(i,j) ) )
                  for( size_t l=0UL; l<N; ++l )
                  {
                     const size_t kbegin( ( IsLower_v<MT2> )
                                          ?( ( IsStrictlyLower_v<MT2> ? l+1UL : l ) )
                                          :( 0UL ) );
                     const size_t kend( ( IsUpper_v<MT2> )
                                        ?( IsStrictlyUpper_v<MT2> ? l : l+1UL )
                                        :( M ) );
                     BLAZE_INTERNAL_ASSERT( kbegin <= kend, "Invalid loop indices detected" );

                     for( size_t k=kbegin; k<kend; ++k )
                        (*lhs)(i*M+k,j*N+l) += A(i,j) * B(k,l);
                  }
      }
      else
      {
         for( size_t i=0UL; i<A.rows(); ++i )
            for( size_t k=0UL; k<M; ++k )
            {
               const size_t lbegin( ( IsUpper_v<MT2> )
                                    ?( ( IsStrictlyUpper_v<MT2> ? k+1UL : k ) )
                                    :( 0UL ) );
               const size_t lend( ( IsLower_v<MT2> )
                                  ?( IsStrictlyLower_v<MT2> ? k : k+1UL )
                                  :( N ) );
               BLAZE_INTERNAL_ASSERT( lbegin <= lend, "Invalid loop indices detected" );

               for( size_t j=0UL; j<A.columns(); ++j )
                  if( !isDefault<strict>( A(i,j) ) )
                     for( size_t l=lbegin; l<lend; ++l )
                        (*lhs)(i*M+k,j*N+l) += A(i,j) * B(k,l);
            }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to column-major dense matrices******************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a dense matrix-dense matrix Kronecker product to a column-major
   //        dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Kronecker product expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense matrix-
   // dense matrix Kronecker product expression to a column-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline void addAssign( DenseMatrix<MT,true>& lhs, const DMatDMatKronExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( rhs.rows() == 0UL || rhs.columns() == 0UL ) {
         return;
      }

      CT1 A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense matrix operand
      CT2 B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side dense matrix operand

      const size_t M( B.rows()    );
      const size_t N( B.columns() );

      if( SO )
      {
         for( size_t j=0UL; j<A.columns(); ++j )
            for( size_t l=0UL; l<N; ++l )
            {
               const size_t kbegin( ( IsLower_v<MT2> )
                                    ?( ( IsStrictlyLower_v<MT2> ? l+1UL : l ) )
                                    :( 0UL ) );
               const size_t kend( ( IsUpper_v<MT2> )
                                  ?( IsStrictlyUpper_v<MT2> ? l : l+1UL )
                                  :( M ) );
               BLAZE_INTERNAL_ASSERT( kbegin <= kend, "Invalid loop indices detected" );

               for( size_t i=0UL; i<A.rows(); ++i )
                  if( !isDefault<strict>( A(i,j) ) )
                     for( size_t k=kbegin; k<kend; ++k )
                        (*lhs)(i*M+k,j*N+l) += A(i,j) * B(k,l);
            }
      }
      else
      {
         for( size_t j=0UL; j<A.columns(); ++j )
            for( size_t i=0UL; i<A.rows(); ++i )
               if( !isDefault<strict>( A(i,j) ) )
                  for( size_t k=0UL; k<M; ++k )
                  {
                     const size_t lbegin( ( IsUpper_v<MT2> )
                                          ?( ( IsStrictlyUpper_v<MT2> ? k+1UL : k ) )
                                          :( 0UL ) );
                     const size_t lend( ( IsLower_v<MT2> )
                                        ?( IsStrictlyLower_v<MT2> ? k : k+1UL )
                                        :( N ) );
                     BLAZE_INTERNAL_ASSERT( lbegin <= lend, "Invalid loop indices detected" );

                     for( size_t l=lbegin; l<lend; ++l )
                        (*lhs)(i*M+k,j*N+l) += A(i,j) * B(k,l);
                  }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse matrices******************************************************
   // No special implementation for the addition assignment to sparse matrices.
   //**********************************************************************************************

   //**Subtraction assignment to row-major dense matrices******************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a dense matrix-dense matrix Kronecker product to a row-major
   //        dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Kronecker product expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtracted assignment of a row-major
   // dense matrix-dense matrix Kronecker product expression to a dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline void subAssign( DenseMatrix<MT,false>& lhs, const DMatDMatKronExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( rhs.rows() == 0UL || rhs.columns() == 0UL ) {
         return;
      }

      CT1 A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense matrix operand
      CT2 B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side dense matrix operand

      const size_t M( B.rows()    );
      const size_t N( B.columns() );

      if( SO )
      {
         for( size_t i=0UL; i<A.rows(); ++i )
            for( size_t j=0UL; j<A.columns(); ++j )
               if( !isDefault<strict>( A(i,j) ) )
                  for( size_t l=0UL; l<N; ++l )
                  {
                     const size_t kbegin( ( IsLower_v<MT2> )
                                          ?( ( IsStrictlyLower_v<MT2> ? l+1UL : l ) )
                                          :( 0UL ) );
                     const size_t kend( ( IsUpper_v<MT2> )
                                        ?( IsStrictlyUpper_v<MT2> ? l : l+1UL )
                                        :( M ) );
                     BLAZE_INTERNAL_ASSERT( kbegin <= kend, "Invalid loop indices detected" );

                     for( size_t k=kbegin; k<kend; ++k )
                        (*lhs)(i*M+k,j*N+l) -= A(i,j) * B(k,l);
                  }
      }
      else
      {
         for( size_t i=0UL; i<A.rows(); ++i )
            for( size_t k=0UL; k<M; ++k )
            {
               const size_t lbegin( ( IsUpper_v<MT2> )
                                    ?( ( IsStrictlyUpper_v<MT2> ? k+1UL : k ) )
                                    :( 0UL ) );
               const size_t lend( ( IsLower_v<MT2> )
                                  ?( IsStrictlyLower_v<MT2> ? k : k+1UL )
                                  :( N ) );
               BLAZE_INTERNAL_ASSERT( lbegin <= lend, "Invalid loop indices detected" );

               for( size_t j=0UL; j<A.columns(); ++j )
                  if( !isDefault<strict>( A(i,j) ) )
                     for( size_t l=lbegin; l<lend; ++l )
                        (*lhs)(i*M+k,j*N+l) -= A(i,j) * B(k,l);
            }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to column-major dense matrices***************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a dense matrix-dense matrix Kronecker product to a
   //        column-major dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Kronecker product expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtracted assignment of a column-major
   // dense matrix-dense matrix Kronecker product expression to a dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline void subAssign( DenseMatrix<MT,true>& lhs, const DMatDMatKronExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( rhs.rows() == 0UL || rhs.columns() == 0UL ) {
         return;
      }

      CT1 A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense matrix operand
      CT2 B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side dense matrix operand

      const size_t M( B.rows()    );
      const size_t N( B.columns() );

      if( SO )
      {
         for( size_t j=0UL; j<A.columns(); ++j )
            for( size_t l=0UL; l<N; ++l )
            {
               const size_t kbegin( ( IsLower_v<MT2> )
                                    ?( ( IsStrictlyLower_v<MT2> ? l+1UL : l ) )
                                    :( 0UL ) );
               const size_t kend( ( IsUpper_v<MT2> )
                                  ?( IsStrictlyUpper_v<MT2> ? l : l+1UL )
                                  :( M ) );
               BLAZE_INTERNAL_ASSERT( kbegin <= kend, "Invalid loop indices detected" );

               for( size_t i=0UL; i<A.rows(); ++i )
                  if( !isDefault<strict>( A(i,j) ) )
                     for( size_t k=kbegin; k<kend; ++k )
                        (*lhs)(i*M+k,j*N+l) -= A(i,j) * B(k,l);
            }
      }
      else
      {
         for( size_t j=0UL; j<A.columns(); ++j )
            for( size_t i=0UL; i<A.rows(); ++i )
               if( !isDefault<strict>( A(i,j) ) )
                  for( size_t k=0UL; k<M; ++k )
                  {
                     const size_t lbegin( ( IsUpper_v<MT2> )
                                          ?( ( IsStrictlyUpper_v<MT2> ? k+1UL : k ) )
                                          :( 0UL ) );
                     const size_t lend( ( IsLower_v<MT2> )
                                        ?( IsStrictlyLower_v<MT2> ? k : k+1UL )
                                        :( N ) );
                     BLAZE_INTERNAL_ASSERT( lbegin <= lend, "Invalid loop indices detected" );

                     for( size_t l=lbegin; l<lend; ++l )
                        (*lhs)(i*M+k,j*N+l) -= A(i,j) * B(k,l);
                  }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse matrices***************************************************
   // No special implementation for the subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**Schur product assignment to row-major dense matrices****************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Schur product assignment of a dense matrix-dense matrix Kronecker product to a
   //        row-major dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Kronecker product expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a dense
   // matrix-dense matrix Kronecker product expression to a row-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline void schurAssign( DenseMatrix<MT,false>& lhs, const DMatDMatKronExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( rhs.rows() == 0UL || rhs.columns() == 0UL ) {
         return;
      }

      CT1 A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense matrix operand
      CT2 B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side dense matrix operand

      const size_t M( B.rows()    );
      const size_t N( B.columns() );

      if( SO )
      {
         for( size_t i=0UL; i<A.rows(); ++i ) {
            for( size_t j=0UL; j<A.columns(); ++j ) {
               if( !isDefault<strict>( A(i,j) ) ) {
                  for( size_t l=0UL; l<N; ++l )
                     for( size_t k=0UL; k<M; ++k )
                        (*lhs)(i*M+k,j*N+l) *= A(i,j) * B(k,l);
               }
               else {
                  for( size_t l=0UL; l<N; ++l )
                     for( size_t k=0UL; k<M; ++k )
                        reset( (*lhs)(i*M+k,j*N+l) );
               }
            }
         }
      }
      else
      {
         for( size_t i=0UL; i<A.rows(); ++i ) {
            for( size_t k=0UL; k<M; ++k ) {
               for( size_t j=0UL; j<A.columns(); ++j ) {
                  if( !isDefault<strict>( A(i,j) ) ) {
                     for( size_t l=0UL; l<N; ++l )
                        (*lhs)(i*M+k,j*N+l) *= A(i,j) * B(k,l);
                  }
                  else {
                     for( size_t l=0UL; l<N; ++l )
                        reset( (*lhs)(i*M+k,j*N+l) );
                  }
               }
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Schur product assignment to column-major dense matrices*************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Schur product assignment of a dense matrix-dense matrix Kronecker product to a
   //        column-major dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Kronecker product expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a dense
   // matrix-dense matrix Kronecker product expression to a column-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline void schurAssign( DenseMatrix<MT,true>& lhs, const DMatDMatKronExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( rhs.rows() == 0UL || rhs.columns() == 0UL ) {
         return;
      }

      CT1 A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense matrix operand
      CT2 B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side dense matrix operand

      const size_t M( B.rows()    );
      const size_t N( B.columns() );

      if( SO )
      {
         for( size_t j=0UL; j<A.columns(); ++j ) {
            for( size_t l=0UL; l<N; ++l ) {
               for( size_t i=0UL; i<A.rows(); ++i ) {
                  if( !isDefault<strict>( A(i,j) ) ) {
                     for( size_t k=0UL; k<M; ++k )
                        (*lhs)(i*M+k,j*N+l) *= A(i,j) * B(k,l);
                  }
                  else {
                     for( size_t k=0UL; k<M; ++k )
                        reset( (*lhs)(i*M+k,j*N+l) );
                  }
               }
            }
         }
      }
      else
      {
         for( size_t j=0UL; j<A.columns(); ++j ) {
            for( size_t i=0UL; i<A.rows(); ++i ) {
               if( !isDefault<strict>( A(i,j) ) ) {
                  for( size_t k=0UL; k<M; ++k )
                     for( size_t l=0UL; l<N; ++l )
                        (*lhs)(i*M+k,j*N+l) *= A(i,j) * B(k,l);
               }
               else {
                  for( size_t k=0UL; k<M; ++k )
                     for( size_t l=0UL; l<N; ++l )
                        reset( (*lhs)(i*M+k,j*N+l) );
               }
            }
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
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( MT2, SO );
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
/*!\brief Computes the Kronecker product of two dense matrices (\f$ A=B \otimes C \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the Kronecker product.
// \param rhs The right-hand side dense matrix for the Kronecker product.
// \return The Kronecker product of the two matrices.
//
// The kron() function computes the Kronecker product of the two given dense matrices:

   \code
   blaze::DynamicMatrix<double> A, B, C;
   // ... Resizing and initialization
   C = kron( A, B );
   \endcode

// The function returns an expression representing a dense matrix of the higher-order element
// type of the two involved matrix element types \a MT1::ElementType and \a MT2::ElementType.
// Both matrix types \a MT1 and \a MT2 as well as the two element types \a MT1::ElementType
// and \a MT2::ElementType have to be supported by the MultTrait class template.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline decltype(auto)
   kron( const DenseMatrix<MT1,SO1>& lhs, const DenseMatrix<MT2,SO2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const DMatDMatKronExpr<MT1,MT2,SO2>;
   return ReturnType( *lhs, *rhs );
}
//*************************************************************************************************

} // namespace blaze

#endif
