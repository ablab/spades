//=================================================================================================
/*!
//  \file blaze/math/expressions/SVecSVecOuterExpr.h
//  \brief Header file for the sparse vector/sparse vector outer product expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_SVECSVECOUTEREXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_SVECSVECOUTEREXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/ColumnVector.h>
#include <blaze/math/constraints/RowVector.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/constraints/Symmetric.h>
#include <blaze/math/constraints/VecTVecMultExpr.h>
#include <blaze/math/constraints/Zero.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/expressions/VecTVecMultExpr.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsTemporary.h>
#include <blaze/math/typetraits/IsZero.h>
#include <blaze/system/MacroDisable.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS SVECSVECOUTEREXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for sparse vector-sparse vector outer products.
// \ingroup sparse_matrix_expression
//
// The SVecSVecOuterExpr class represents the compile time expression for sparse vector-sparse
// vector outer products.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side sparse vector
class SVecSVecOuterExpr
   : public VecTVecMultExpr< SparseMatrix< SVecSVecOuterExpr<VT1,VT2>, false > >
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   using RT1 = ResultType_t<VT1>;     //!< Result type of the left-hand side sparse vector expression.
   using RT2 = ResultType_t<VT2>;     //!< Result type of the right-hand side sparse vector expression.
   using RN1 = ReturnType_t<VT1>;     //!< Return type of the left-hand side sparse vector expression.
   using RN2 = ReturnType_t<VT2>;     //!< Return type of the right-hand side sparse vector expression.
   using CT1 = CompositeType_t<VT1>;  //!< Composite type of the left-hand side sparse vector expression.
   using CT2 = CompositeType_t<VT2>;  //!< Composite type of the right-hand side sparse vector expression.
   //**********************************************************************************************

   //**Return type evaluation**********************************************************************
   //! Compilation switch for the selection of the subscript operator return type.
   /*! The \a returnExpr compile time constant expression is a compilation switch for the
       selection of the \a ReturnType. If either vector operand returns a temporary vector
       or matrix, \a returnExpr will be set to \a false and the subscript operator will
       return it's result by value. Otherwise \a returnExpr will be set to \a true and
       the subscript operator may return it's result as an expression. */
   static constexpr bool returnExpr = ( !IsTemporary_v<RN1> && !IsTemporary_v<RN2> );

   //! Expression return type for the subscript operator.
   using ExprReturnType = decltype( std::declval<RN1>() * std::declval<RN2>() );
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this SVecSVecOuterExpr instance.
   using This = SVecSVecOuterExpr<VT1,VT2>;

   //! Base type of this SVecSVecOuterExpr instance.
   using BaseType = VecTVecMultExpr< SparseMatrix<This,false> >;

   using ResultType    = MultTrait_t<RT1,RT2>;         //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<ResultType>;   //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<ResultType>;    //!< Resulting element type.

   //! Return type for expression template evaluations.
   using ReturnType = const If_t< returnExpr, ExprReturnType, ElementType >;

   //! Data type for composite expression templates.
   using CompositeType = const ResultType;

   //! Composite type of the left-hand side sparse vector expression.
   using LeftOperand = If_t< IsExpression_v<VT1>, const VT1, const VT1& >;

   //! Composite type of the right-hand side sparse vector expression.
   using RightOperand = If_t< IsExpression_v<VT2>, const VT2, const VT2& >;

   //! Type for the assignment of the left-hand side sparse vector operand.
   using LT = If_t< IsComputation_v<VT1>, const RT1, CT1 >;

   //! Type for the assignment of the right-hand side sparse vector operand.
   using RT = If_t< IsComputation_v<VT2>, const RT2, CT2 >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = false;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the SVecSVecOuterExpr class.
   //
   // \param lhs The left-hand side sparse vector operand of the multiplication expression.
   // \param rhs The right-hand side sparse vector operand of the multiplication expression.
   */
   inline SVecSVecOuterExpr( const VT1& lhs, const VT2& rhs ) noexcept
      : lhs_( lhs )  // Left-hand side sparse vector of the multiplication expression
      , rhs_( rhs )  // Right-hand side sparse vector of the multiplication expression
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
      BLAZE_INTERNAL_ASSERT( i < lhs_.size(), "Invalid row access index"    );
      BLAZE_INTERNAL_ASSERT( j < rhs_.size(), "Invalid column access index" );

      return lhs_[i] * rhs_[j];
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
      if( i >= lhs_.size() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
      }
      if( j >= rhs_.size() ) {
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
      return lhs_.size();
   }
   //**********************************************************************************************

   //**Columns function****************************************************************************
   /*!\brief Returns the current number of columns of the matrix.
   //
   // \return The number of columns of the matrix.
   */
   inline size_t columns() const noexcept {
      return rhs_.size();
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
      return ( isDefault( lhs_[i] ) )?( size_t(0) ):( rhs_.nonZeros(i) );
   }
   //**********************************************************************************************

   //**Left operand access*************************************************************************
   /*!\brief Returns the left-hand side sparse vector operand.
   //
   // \return The left-hand side sparse vector operand.
   */
   inline LeftOperand leftOperand() const noexcept {
      return lhs_;
   }
   //**********************************************************************************************

   //**Right operand access************************************************************************
   /*!\brief Returns the right-hand side sparse vector operand.
   //
   // \return The right-hand side sparse vector operand.
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
   LeftOperand  lhs_;  //!< Left-hand side sparse vector of the multiplication expression.
   RightOperand rhs_;  //!< Right-hand side sparse vector of the multiplication expression.
   //**********************************************************************************************

   //**Assignment to row-major dense matrices******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse vector-sparse vector outer product to a row-major dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side outer product expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse vector-sparse
   // vector outer product expression to a row-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline void assign( DenseMatrix<MT,false>& lhs, const SVecSVecOuterExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT x( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side sparse vector operand
      RT y( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (*lhs).rows()   , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == (*lhs).columns(), "Invalid vector size" );

      const auto lend( x.end() );
      const auto rend( y.end() );

      for( auto lelem=x.begin(); lelem!=lend; ++lelem ) {
         if( !isDefault( lelem->value() ) ) {
            for( auto relem=y.begin(); relem!=rend; ++relem ) {
               (*lhs)(lelem->index(),relem->index()) = lelem->value() * relem->value();
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to column-major dense matrices***************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse vector-sparse vector outer product to a column-major
   //        dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side outer product expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse vector-sparse
   // vector outer product expression to a column-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline void assign( DenseMatrix<MT,true>& lhs, const SVecSVecOuterExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT x( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side sparse vector operand
      RT y( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (*lhs).rows()   , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == (*lhs).columns(), "Invalid vector size" );

      const auto lend( x.end() );
      const auto rend( y.end() );

      for( auto relem=y.begin(); relem!=rend; ++relem ) {
         if( !isDefault( relem->value() ) ) {
            for( auto lelem=x.begin(); lelem!=lend; ++lelem ) {
               (*lhs)(lelem->index(),relem->index()) = lelem->value() * relem->value();
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to row-major sparse matrices*****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse vector-sparse vector outer product to a row-major sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side outer product expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse vector-sparse
   // vector outer product expression to a row-major sparse matrix.
   */
   template< typename MT >  // Type of the target sparse matrix
   friend inline void assign( SparseMatrix<MT,false>& lhs, const SVecSVecOuterExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()     == rhs.rows()    , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns()  == rhs.columns() , "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( (*lhs).capacity() >= rhs.nonZeros(), "Insufficient capacity"     );

      LT x( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side sparse vector operand
      RT y( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (*lhs).rows()   , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == (*lhs).columns(), "Invalid vector size" );

      // Final memory allocation (based on the evaluated operands)
      (*lhs).reserve( x.nonZeros() * y.nonZeros() );

      // Performing the outer product
      const auto lend( x.end() );
      const auto rend( y.end() );
      size_t index( 0UL );

      for( auto lelem=x.begin(); lelem!=lend; ++lelem ) {
         if( !isDefault( lelem->value() ) ) {
            for( ; index < lelem->index(); ++index ) {
               (*lhs).finalize( index );
            }
            for( auto relem=y.begin(); relem!=rend; ++relem ) {
               (*lhs).append( lelem->index(), relem->index(), lelem->value() * relem->value() );
            }
            (*lhs).finalize( index++ );
         }
      }

      for( ; index < x.size(); ++index ) {
         (*lhs).finalize( index );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to column-major sparse matrices**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse vector-sparse vector outer product to a column-major
   //        sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side outer product expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse vector-sparse
   // vector outer product expression to a column-major sparse matrix.
   */
   template< typename MT >  // Type of the target sparse matrix
   friend inline void assign( SparseMatrix<MT,true>& lhs, const SVecSVecOuterExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()     == rhs.rows()    , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns()  == rhs.columns() , "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( (*lhs).capacity() >= rhs.nonZeros(), "Insufficient capacity"     );

      LT x( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side sparse vector operand
      RT y( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (*lhs).rows()   , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == (*lhs).columns(), "Invalid vector size" );

      const auto lend( x.end() );
      const auto rend( y.end() );
      size_t index( 0UL );

      for( auto relem=y.begin(); relem!=rend; ++relem ) {
         if( !isDefault( relem->value() ) ) {
            for( ; index < relem->index(); ++index ) {
               (*lhs).finalize( index );
            }
            for( auto lelem=x.begin(); lelem!=lend; ++lelem ) {
               (*lhs).append( lelem->index(), relem->index(), lelem->value() * relem->value() );
            }
            (*lhs).finalize( index++ );
         }
      }

      for( ; index < y.size(); ++index ) {
         (*lhs).finalize( index );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to row-major dense matrices*********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a sparse vector-sparse vector outer product to a row-major
   //        dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side outer product expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a sparse vector-
   // sparse vector outer product expression to a row-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline void addAssign( DenseMatrix<MT,false>& lhs, const SVecSVecOuterExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT x( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side sparse vector operand
      RT y( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (*lhs).rows()   , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == (*lhs).columns(), "Invalid vector size" );

      const auto lend( x.end() );
      const auto rend( y.end() );

      for( auto lelem=x.begin(); lelem!=lend; ++lelem ) {
         if( !isDefault( lelem->value() ) ) {
            for( auto relem=y.begin(); relem!=rend; ++relem ) {
               (*lhs)(lelem->index(),relem->index()) += lelem->value() * relem->value();
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to column-major dense matrices******************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a sparse vector-sparse vector outer product to a column-major
   //        dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side outer product expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a sparse vector-
   // sparse vector outer product expression to a column-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline void addAssign( DenseMatrix<MT,true>& lhs, const SVecSVecOuterExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT x( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side sparse vector operand
      RT y( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (*lhs).rows()   , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == (*lhs).columns(), "Invalid vector size" );

      const auto lend( x.end() );
      const auto rend( y.end() );

      for( auto relem=y.begin(); relem!=rend; ++relem ) {
         if( !isDefault( relem->value() ) ) {
            for( auto lelem=x.begin(); lelem!=lend; ++lelem ) {
               (*lhs)(lelem->index(),relem->index()) += lelem->value() * relem->value();
            }
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
   /*!\brief Subtraction assignment of a sparse vector-sparse vector outer product to a row-major
   //        dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side outer product expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a sparse vector-
   // sparse vector outer product expression to a row-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline void subAssign( DenseMatrix<MT,false>& lhs, const SVecSVecOuterExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT x( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side sparse vector operand
      RT y( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (*lhs).rows()   , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == (*lhs).columns(), "Invalid vector size" );

      const auto lend( x.end() );
      const auto rend( y.end() );

      for( auto lelem=x.begin(); lelem!=lend; ++lelem ) {
         if( !isDefault( lelem->value() ) ) {
            for( auto relem=y.begin(); relem!=rend; ++relem ) {
               (*lhs)(lelem->index(),relem->index()) -= lelem->value() * relem->value();
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to column-major dense matrices***************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a sparse vector-sparse vector outer product to a column-major
   //        dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side outer product expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a sparse vector-
   // sparse vector outer product expression to a column-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline void subAssign( DenseMatrix<MT,true>& lhs, const SVecSVecOuterExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT x( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side sparse vector operand
      RT y( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (*lhs).rows()   , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == (*lhs).columns(), "Invalid vector size" );

      const auto lend( x.end() );
      const auto rend( y.end() );

      for( auto relem=y.begin(); relem!=rend; ++relem ) {
         if( !isDefault( relem->value() ) ) {
            for( auto lelem=x.begin(); lelem!=lend; ++lelem ) {
               (*lhs)(lelem->index(),relem->index()) -= lelem->value() * relem->value();
            }
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
   /*!\brief Schur product assignment of a sparse vector-sparse vector outer product to a row-major
   //        dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side outer product expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a sparse
   // vector-sparse vector outer product expression to a row-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline void schurAssign( DenseMatrix<MT,false>& lhs, const SVecSVecOuterExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT x( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side sparse vector operand
      RT y( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (*lhs).rows()   , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == (*lhs).columns(), "Invalid vector size" );

      const auto lend( x.end() );
      const auto rend( y.end() );

      size_t i( 0UL );

      for( auto lelem=x.begin(); lelem!=lend; ++lelem )
      {
         if( isDefault( lelem->value() ) ) continue;

         for( ; i<lelem->index(); ++i ) {
            for( size_t j=0UL; j<y.size(); ++j )
               reset( (*lhs)(i,j) );
         }

         size_t j( 0UL );

         for( auto relem=y.begin(); relem!=rend; ++relem, ++j ) {
            for( ; j<relem->index(); ++j )
               reset( (*lhs)(i,j) );
            (*lhs)(lelem->index(),relem->index()) *= lelem->value() * relem->value();
         }

         for( ; j<y.size(); ++j ) {
            reset( (*lhs)(i,j) );
         }

         ++i;
      }

      for( ; i<x.size(); ++i ) {
         for( size_t j=0UL; j<y.size(); ++j )
            reset( (*lhs)(i,j) );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Schur product assignment to column-major dense matrices*************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Schur product assignment of a sparse vector-sparse vector outer product to a
   //        column-major dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side outer product expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a sparse
   // vector-sparse vector outer product expression to a column-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline void schurAssign( DenseMatrix<MT,true>& lhs, const SVecSVecOuterExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT x( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side sparse vector operand
      RT y( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (*lhs).rows()   , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == (*lhs).columns(), "Invalid vector size" );

      const auto lend( x.end() );
      const auto rend( y.end() );

      size_t j( 0UL );

      for( auto relem=y.begin(); relem!=rend; ++relem )
      {
         if( isDefault( relem->value() ) ) continue;

         for( ; j<relem->index(); ++j ) {
            for( size_t i=0UL; i<x.size(); ++i )
               reset( (*lhs)(i,j) );
         }

         size_t i( 0UL );

         for( auto lelem=x.begin(); lelem!=lend; ++lelem, ++i ) {
            for( ; i<lelem->index(); ++i )
               reset( (*lhs)(i,j) );
            (*lhs)(lelem->index(),relem->index()) *= lelem->value() * relem->value();
         }

         for( ; i<x.size(); ++i ) {
            reset( (*lhs)(i,j) );
         }

         ++j;
      }

      for( ; j<y.size(); ++j ) {
         for( size_t i=0UL; i<x.size(); ++i )
            reset( (*lhs)(i,j) );
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
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( VT1 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( VT1 );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( VT2 );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE   ( VT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE     ( VT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE     ( VT2 );
   BLAZE_CONSTRAINT_MUST_FORM_VALID_VECTVECMULTEXPR( VT1, VT2 );
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
/*!\brief Backend implementation of the sparse vector-sparse vector outer product
//        (\f$ A=\vec{b}*\vec{c}^T \f$).
// \ingroup sparse_matrix
//
// \param lhs The left-hand side sparse vector for the outer product.
// \param rhs The right-hand side transpose sparse vector for the outer product.
// \return The resulting sparse matrix.
//
// This function implements a performance optimized treatment of the sparse vector-sparse vector
// outer product.
*/
template< typename VT1  // Type of the left-hand side sparse vector
        , typename VT2  // Type of the right-hand side sparse vector
        , DisableIf_t< IsZero_v<VT1> || IsZero_v<VT2> >* = nullptr >
inline const SVecSVecOuterExpr<VT1,VT2>
   svecsvecouter( const SparseVector<VT1,false>& lhs, const SparseVector<VT2,true>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return SVecSVecOuterExpr<VT1,VT2>( *lhs, *rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the (zero) sparse vector-(zero) sparse vector outer product
//        (\f$ A=\vec{b}*\vec{c}^T \f$).
// \ingroup sparse_matrix
//
// \param lhs The left-hand side sparse vector for the outer product.
// \param rhs The right-hand side transpose sparse vector for the outer product.
// \return The zero sparse matrix.
//
// This function implements a performance optimized treatment of the (zero) sparse vector-(zero)
// sparse vector outer product. It returns a zero matrix.
*/
template< typename VT1  // Type of the left-hand side sparse vector
        , typename VT2  // Type of the right-hand side sparse vector
        , EnableIf_t< IsZero_v<VT1> || IsZero_v<VT2> >* = nullptr >
inline decltype(auto)
   svecsvecouter( const SparseVector<VT1,false>& lhs, const SparseVector<VT2,true>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const MultTrait_t< ResultType_t<VT1>, ResultType_t<VT2> >;

   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( ReturnType );
   BLAZE_CONSTRAINT_MUST_BE_ZERO_TYPE( ReturnType );

   return ReturnType( (*lhs).size(), (*rhs).size() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the sparse vector-sparse vector outer product
//        (\f$ A=\vec{b}*\vec{c}^T \f$).
// \ingroup sparse_matrix
//
// \param lhs The left-hand side sparse vector for the outer product.
// \param rhs The right-hand side transpose sparse vector for the outer product.
// \return The resulting sparse matrix.
//
// This operator represents the outer product between a sparse vector and a transpose sparse
// vector:

   \code
   using blaze::columnVector;
   using blaze::rowMajor;

   blaze::CompressedVector<double,columnVector> a, b;
   blaze::CompressedMatrix<double,rowMajor> A;
   // ... Resizing and initialization
   A = a * trans(b);
   \endcode

// The operator returns an expression representing a sparse matrix of the higher-order element
// type of the two involved element types \a VT1::ElementType and \a VT2::ElementType. Both
// vector types \a VT1 and \a VT2 as well as the two element types \a VT1::ElementType and
// \a VT2::ElementType have to be supported by the MultTrait class template.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side sparse vector
inline decltype(auto)
   operator*( const SparseVector<VT1,false>& lhs, const SparseVector<VT2,true>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return svecsvecouter( *lhs, *rhs );
}
//*************************************************************************************************

} // namespace blaze

#endif
