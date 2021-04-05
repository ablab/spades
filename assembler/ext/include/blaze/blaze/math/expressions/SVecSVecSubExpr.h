//=================================================================================================
/*!
//  \file blaze/math/expressions/SVecSVecSubExpr.h
//  \brief Header file for the sparse vector/sparse vector subtraction expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_SVECSVECSUBEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_SVECSVECSUBEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/constraints/VecVecSubExpr.h>
#include <blaze/math/constraints/Zero.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/expressions/VecVecSubExpr.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsTemporary.h>
#include <blaze/math/typetraits/IsZero.h>
#include <blaze/system/MacroDisable.h>
#include <blaze/util/algorithms/Min.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsSame.h>


namespace blaze {

//=================================================================================================
//
//  CLASS SVECSVECSUBEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for sparse vector-sparse vector subtractions.
// \ingroup sparse_vector_expression
//
// The SVecSVecSubExpr class represents the compile time expression for subtractions between
// sparse vectors.
*/
template< typename VT1  // Type of the left-hand side sparse vector
        , typename VT2  // Type of the right-hand side sparse vector
        , bool TF >     // Transpose flag
class SVecSVecSubExpr
   : public VecVecSubExpr< SparseVector< SVecSVecSubExpr<VT1,VT2,TF>, TF > >
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
   using ExprReturnType = decltype( std::declval<RN1>() - std::declval<RN2>() );
   //**********************************************************************************************

   //**Parallel evaluation strategy****************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! This variable template is a helper for the selection of the parallel evaluation strategy.
       In case the target vector is SMP assignable, the variable is set to 1 and the expression
       specific evaluation strategy is selected. Otherwise the variable is set to 0 and the
       default strategy is chosen. */
   template< typename VT >
   static constexpr bool UseSMPAssign_v = VT::smpAssignable;
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this SVecSVecSubExpr instance.
   using This = SVecSVecSubExpr<VT1,VT2,TF>;

   //! Base type of this SVecSVecSubExpr instance.
   using BaseType = VecVecSubExpr< SparseVector<This,TF> >;

   using ResultType    = SubTrait_t<RT1,RT2>;          //!< Result type for expression template evaluations.
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
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = false;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the SVecSVecSubExpr class.
   */
   inline SVecSVecSubExpr( const VT1& lhs, const VT2& rhs ) noexcept
      : lhs_( lhs )  // Left-hand side sparse vector of the subtraction expression
      , rhs_( rhs )  // Right-hand side sparse vector of the subtraction expression
   {
      BLAZE_INTERNAL_ASSERT( lhs.size() == rhs.size(), "Invalid vector sizes" );
   }
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < lhs_.size(), "Invalid vector access index" );
      return lhs_[index] - rhs_[index];
   }
   //**********************************************************************************************

   //**At function*********************************************************************************
   /*!\brief Checked access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   // \exception std::out_of_range Invalid vector access index.
   */
   inline ReturnType at( size_t index ) const {
      if( index >= lhs_.size() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid vector access index" );
      }
      return (*this)[index];
   }
   //**********************************************************************************************

   //**Size function*******************************************************************************
   /*!\brief Returns the current size/dimension of the vector.
   //
   // \return The size of the vector.
   */
   inline size_t size() const noexcept {
      return lhs_.size();
   }
   //**********************************************************************************************

   //**NonZeros function***************************************************************************
   /*!\brief Returns the number of non-zero elements in the sparse vector.
   //
   // \return The number of non-zero elements in the sparse vector.
   */
   inline size_t nonZeros() const {
      return min( lhs_.size(), lhs_.nonZeros() + rhs_.nonZeros() );
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
   LeftOperand  lhs_;  //!< Left-hand side sparse vector of the subtraction expression.
   RightOperand rhs_;  //!< Right-hand side sparse vector of the subtraction expression.
   //**********************************************************************************************

   //**Default assignment to dense vectors*********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default assignment of a sparse vector-sparse vector subtraction to a dense vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side subtraction expression to be assigned.
   // \return void
   //
   // This function implements the default assignment of a sparse vector-sparse vector
   // subtraction expression to a dense vector. This function is used in case the element
   // type is resizable.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline auto assign( DenseVector<VT,TF>& lhs, const SVecSVecSubExpr& rhs )
      -> EnableIf_t< IsResizable_v< ElementType_t<VT> > >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      CT1 x( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side sparse vector operand
      CT2 y( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size(), "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size(), "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (*lhs).size()  , "Invalid vector size" );

      const auto lend( x.end() );
      const auto rend( y.end() );

      for( auto l=x.begin(); l!=lend; ++l ) {
         (*lhs)[l->index()] = l->value();
      }

      for( auto r=y.begin(); r!=rend; ++r ) {
         if( isDefault( (*lhs)[r->index()] ) )
            (*lhs)[r->index()] = -r->value();
         else
            (*lhs)[r->index()] -= r->value();
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Optimized assignment to dense vectors*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Optimized assignment of a sparse vector-sparse vector subtraction to a dense vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side subtraction expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse vector-sparse
   // vector subtraction expression to a dense vector. This function is used in case the element
   // type is not resizable.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline auto assign( DenseVector<VT,TF>& lhs, const SVecSVecSubExpr& rhs )
      -> DisableIf_t< IsResizable_v< ElementType_t<VT> > >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      CT1 x( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side sparse vector operand
      CT2 y( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size(), "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size(), "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (*lhs).size()  , "Invalid vector size" );

      const auto lend( x.end() );
      const auto rend( y.end() );

      for( auto l=x.begin(); l!=lend; ++l ) {
         (*lhs)[l->index()] = l->value();
      }

      for( auto r=y.begin(); r!=rend; ++r ) {
         (*lhs)[r->index()] -= r->value();
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse vectors****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse vector-sparse vector subtraction to a sparse vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side subtraction expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse vector-sparse
   // vector subtraction expression to a sparse vector.
   */
   template< typename VT >  // Type of the target sparse vector
   friend inline void assign( SparseVector<VT,TF>& lhs, const SVecSVecSubExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      CT1 x( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side sparse vector operand
      CT2 y( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size(), "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size(), "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (*lhs).size()  , "Invalid vector size" );

      // Final memory allocation (based on the evaluated operands)
      (*lhs).reserve( min( x.size(), x.nonZeros() + y.nonZeros() ) );

      // Performing the vector subtraction
      const auto lend( x.end() );
      const auto rend( y.end() );

      auto l( x.begin() );
      auto r( y.begin() );

      while( l != lend && r != rend )
      {
         if( l->index() < r->index() ) {
            (*lhs).append( l->index(), l->value() );
            ++l;
         }
         else if( l->index() > r->index() ) {
            (*lhs).append( r->index(), -r->value() );
            ++r;
         }
         else {
            (*lhs).append( l->index(), l->value() - r->value() );
            ++l;
            ++r;
         }
      }

      while( l != lend ) {
         (*lhs).append( l->index(), l->value() );
         ++l;
      }

      while( r != rend ) {
         (*lhs).append( r->index(), -r->value() );
         ++r;
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense vectors********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a sparse vector-sparse vector subtraction to a dense vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side subtraction expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a sparse vector-
   // sparse vector subtraction expression to a dense vector.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline void addAssign( DenseVector<VT,TF>& lhs, const SVecSVecSubExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      addAssign( *lhs, rhs.lhs_ );
      subAssign( *lhs, rhs.rhs_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse vectors*******************************************************
   // No special implementation for the addition assignment to sparse vectors.
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors*****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a sparse vector-sparse vector subtraction to a dense vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side subtraction expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a sparse vector-
   // sparse vector subtraction expression to a dense vector.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline void subAssign( DenseVector<VT,TF>& lhs, const SVecSVecSubExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      subAssign( *lhs, rhs.lhs_ );
      addAssign( *lhs, rhs.rhs_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse vectors****************************************************
   // No special implementation for the subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**Multiplication assignment to dense vectors**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Multiplication assignment of a sparse vector-sparse vector subtraction to a dense vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side subtraction expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a sparse
   // vector-sparse vector subtraction expression to a dense vector.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline void multAssign( DenseVector<VT,TF>& lhs, const SVecSVecSubExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      multAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Multiplication assignment to sparse vectors*************************************************
   // No special implementation for the multiplication assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP assignment to dense vectors*************************************************************
   // No special implementation for the SMP assignment to dense vectors.
   //**********************************************************************************************

   //**SMP assignment to sparse vectors************************************************************
   // No special implementation for the SMP assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP addition assignment to dense vectors****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a sparse vector-sparse vector subtraction to a dense vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side subtraction expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a sparse
   // vector-sparse vector subtraction expression to a dense vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline auto smpAddAssign( DenseVector<VT,TF>& lhs, const SVecSVecSubExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      smpAddAssign( *lhs, rhs.lhs_ );
      smpSubAssign( *lhs, rhs.rhs_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to sparse vectors***************************************************
   // No special implementation for the SMP addition assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP subtraction assignment to dense vectors*************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a sparse vector-sparse vector subtraction to a dense
   //        vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side subtraction expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a sparse
   // vector-sparse vector subtraction expression to a dense vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline auto smpSubAssign( DenseVector<VT,TF>& lhs, const SVecSVecSubExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      smpSubAssign( *lhs, rhs.lhs_ );
      smpAddAssign( *lhs, rhs.rhs_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to sparse vectors************************************************
   // No special implementation for the SMP subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP multiplication assignment to dense vectors**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP multiplication assignment of a sparse vector-sparse vector subtraction to a dense
   //        vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side subtraction expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized SMP multiplication assignment of a sparse
   // vector-sparse vector subtraction expression to a dense vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline auto smpMultAssign( DenseVector<VT,TF>& lhs, const SVecSVecSubExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( rhs );
      smpMultAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP multiplication assignment to sparse vectors*********************************************
   // No special implementation for the SMP multiplication assignment to sparse vectors.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( VT1 );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( VT2 );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT1, TF );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT2, TF );
   BLAZE_CONSTRAINT_MUST_FORM_VALID_VECVECSUBEXPR( VT1, VT2 );
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
/*!\brief Backend implementation of the subtraction of two sparse vectors (\f$ \vec{a}=\vec{b}-\vec{c} \f$).
// \ingroup sparse_vector
//
// \param lhs The left-hand side sparse vector for the subtraction.
// \param rhs The right-hand side sparse vector for the subtraction.
// \return The difference of the two vectors.
//
// This function implements a performance optimized treatment of the subtraction between two
// sparse vectors.
*/
template< typename VT1  // Type of the left-hand side sparse vector
        , typename VT2  // Type of the right-hand side sparse vector
        , bool TF       // Transpose flag
        , DisableIf_t< ( ( IsZero_v<VT1> || IsZero_v<VT2> ) &&
                         IsSame_v< ElementType_t<VT1>, ElementType_t<VT2> > ) ||
                       ( IsZero_v<VT1> && IsZero_v<VT2> ) >* = nullptr >
inline const SVecSVecSubExpr<VT1,VT2,TF>
   svecsvecsub( const SparseVector<VT1,TF>& lhs, const SparseVector<VT2,TF>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).size() == (*rhs).size(), "Invalid vector sizes" );

   return SVecSVecSubExpr<VT1,VT2,TF>( *lhs, *rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the subtraction between a sparse vector and a zero vector
//        (\f$ \vec{a}=\vec{b}-\vec{c} \f$).
// \ingroup sparse_vector
//
// \param lhs The left-hand side sparse vector for the subtraction.
// \param rhs The right-hand side zero vector for the subtraction.
// \return The difference of the two vectors.
//
// This function implements a performance optimized treatment of the subtraction between a
// sparse vector and a zero vector.
*/
template< typename VT1  // Type of the left-hand side sparse vector
        , typename VT2  // Type of the right-hand side sparse vector
        , bool TF       // Transpose flag
        , EnableIf_t< !IsZero_v<VT1> && IsZero_v<VT2> &&
                      IsSame_v< ElementType_t<VT1>, ElementType_t<VT2> > >* = nullptr >
inline const VT1&
   svecsvecsub( const SparseVector<VT1,TF>& lhs, const SparseVector<VT2,TF>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( rhs );

   BLAZE_INTERNAL_ASSERT( (*lhs).size() == (*rhs).size(), "Invalid vector sizes" );

   return (*lhs);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the subtraction between a zero vector and a sparse vector
//        (\f$ \vec{a}=\vec{b}-\vec{c} \f$).
// \ingroup sparse_vector
//
// \param lhs The left-hand side zero vector for the subtraction.
// \param rhs The right-hand side sparse vector for the subtraction.
// \return The difference of the two vectors.
//
// This function implements a performance optimized treatment of the subtraction between a zero
// vector and a sparse vector.
*/
template< typename VT1  // Type of the left-hand side sparse vector
        , typename VT2  // Type of the right-hand side sparse vector
        , bool TF       // Transpose flag
        , EnableIf_t< IsZero_v<VT1> && !IsZero_v<VT2> &&
                      IsSame_v< ElementType_t<VT1>, ElementType_t<VT2> > >* = nullptr >
inline decltype(auto)
   svecsvecsub( const SparseVector<VT1,TF>& lhs, const SparseVector<VT2,TF>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( lhs );

   BLAZE_INTERNAL_ASSERT( (*lhs).size() == (*rhs).size(), "Invalid vector sizes" );

   return -(*rhs);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the subtraction of two zero vectors (\f$ \vec{a}=\vec{b}-\vec{c} \f$).
// \ingroup sparse_vector
//
// \param lhs The left-hand side zero vector for the subtraction.
// \param rhs The right-hand side zero vector for the subtraction.
// \return The resulting zero vector.
//
// This function implements a performance optimized treatment of the subtraction between two zero
// vectors. It returns a zero vector.
*/
template< typename VT1  // Type of the left-hand side sparse vector
        , typename VT2  // Type of the right-hand side sparse vector
        , bool TF       // Transpose flag
        , EnableIf_t< IsZero_v<VT1> && IsZero_v<VT2> >* = nullptr >
inline decltype(auto)
   svecsvecsub( const SparseVector<VT1,TF>& lhs, const SparseVector<VT2,TF>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( rhs );

   BLAZE_INTERNAL_ASSERT( (*lhs).size() == (*rhs).size(), "Invalid vector sizes" );

   using ReturnType = const SubTrait_t< ResultType_t<VT1>, ResultType_t<VT2> >;

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ReturnType, TF );
   BLAZE_CONSTRAINT_MUST_BE_ZERO_TYPE( ReturnType );

   return ReturnType( (*lhs).size() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction operator for the subtraction of two sparse vectors (\f$ \vec{a}=\vec{b}-\vec{c} \f$).
// \ingroup sparse_vector
//
// \param lhs The left-hand side sparse vector for the vector subtraction.
// \param rhs The right-hand side sparse vector to be subtracted from the vector.
// \return The difference of the two sparse vectors.
// \exception std::invalid_argument Vector sizes do not match.
//
// This operator represents the subtraction of two sparse vectors:

   \code
   blaze::CompressedVector<double> a, b, c;
   // ... Resizing and initialization
   c = a - b;
   \endcode

// The operator returns a sparse vector of the higher-order element type of the two involved
// vector element types \a VT1::ElementType and \a VT2::ElementType. Both vector types \a VT1
// and \a VT2 as well as the two element types \a VT1::ElementType and \a VT2::ElementType
// have to be supported by the SubTrait class template.\n
// In case the current sizes of the two given vectors don't match, a \a std::invalid_argument
// is thrown.
*/
template< typename VT1  // Type of the left-hand side sparse vector
        , typename VT2  // Type of the right-hand side sparse vector
        , bool TF >     // Transpose flag
inline decltype(auto)
   operator-( const SparseVector<VT1,TF>& lhs, const SparseVector<VT2,TF>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   if( (*lhs).size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   return svecsvecsub( *lhs, *rhs );
}
//*************************************************************************************************

} // namespace blaze

#endif
