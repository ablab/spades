//=================================================================================================
/*!
//  \file blaze/math/expressions/DVecSVecAddExpr.h
//  \brief Header file for the dense vector/sparse vector addition expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DVECSVECADDEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DVECSVECADDEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/constraints/VecVecAddExpr.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/VecVecAddExpr.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsTemporary.h>
#include <blaze/math/typetraits/IsZero.h>
#include <blaze/system/MacroDisable.h>
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
//  CLASS DVECSVECADDEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for dense vector-sparse vector additions.
// \ingroup dense_vector_expression
//
// The DVecSVecAddExpr class represents the compile time expression for additions between
// a dense vector and a sparse vector.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side sparse vector
        , bool TF >     // Transpose flag
class DVecSVecAddExpr
   : public VecVecAddExpr< DenseVector< DVecSVecAddExpr<VT1,VT2,TF>, TF > >
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   using RT1 = ResultType_t<VT1>;     //!< Result type of the left-hand side dense vector expression.
   using RT2 = ResultType_t<VT2>;     //!< Result type of the right-hand side sparse vector expression.
   using RN1 = ReturnType_t<VT1>;     //!< Return type of the left-hand side dense vector expression.
   using RN2 = ReturnType_t<VT2>;     //!< Return type of the right-hand side sparse vector expression.
   using CT1 = CompositeType_t<VT1>;  //!< Composite type of the left-hand side dense vector expression.
   using CT2 = CompositeType_t<VT2>;  //!< Composite type of the right-hand side sparse vector expression.
   using TT1 = TransposeType_t<VT1>;  //!< Transpose type of the left-hand side dense vector expression.
   using TT2 = TransposeType_t<VT2>;  //!< Transpose type of the right-hand side sparse vector expression.
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
   using ExprReturnType = decltype( std::declval<RN1>() + std::declval<RN2>() );
   //**********************************************************************************************

   //**Parallel evaluation strategy****************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! This variable template is a helper for the selection of the parallel evaluation strategy.
       In case at least one of the two vector operands is not SMP assignable, the variable is set
       to 1 and the expression specific evaluation strategy is selected. Otherwise the variable
       is set to 0 and the default strategy is chosen. */
   template< typename VT >
   static constexpr bool UseSMPAssign_v = ( !VT1::smpAssignable || !VT2::smpAssignable );
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this DVecSVecAddExpr instance.
   using This = DVecSVecAddExpr<VT1,VT2,TF>;

   //! Base type of this DVecSVecAddExpr instance.
   using BaseType = VecVecAddExpr< DenseVector<This,TF> >;

   using ResultType    = AddTrait_t<RT1,RT2>;          //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<ResultType>;    //!< Resulting element type.

   //! Return type for expression template evaluations.
   using ReturnType = const If_t< returnExpr, ExprReturnType, ElementType >;

   //! Data type for composite expression templates.
   using CompositeType = const ResultType;

   //! Composite type of the left-hand side dense vector expression.
   using LeftOperand = If_t< IsExpression_v<VT1>, const VT1, const VT1& >;

   //! Composite type of the right-hand side sparse vector expression.
   using RightOperand = If_t< IsExpression_v<VT2>, const VT2, const VT2& >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = false;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = false;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DVecSVecAddExpr class.
   //
   // \param lhs The left-hand side dense vector operand of the addition expression.
   // \param rhs The right-hand side sparse vector operand of the addition expression.
   */
   inline DVecSVecAddExpr( const VT1& lhs, const VT2& rhs ) noexcept
      : lhs_( lhs )  // Left-hand side dense vector of the addition expression
      , rhs_( rhs )  // Right-hand side sparse vector of the addition expression
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
      return lhs_[index] + rhs_[index];
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

   //**Left operand access*************************************************************************
   /*!\brief Returns the left-hand side dense vector operand.
   //
   // \return The left-hand side dense vector operand.
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
   // \return \a true in case the expression can alias effect is detected, \a false otherwise.
   */
   template< typename T >
   inline bool canAlias( const T* alias ) const noexcept {
      return ( IsExpression_v<VT1> && lhs_.canAlias( alias ) ) ||
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
      return ( lhs_.isAliased( alias ) || rhs_.isAliased( alias ) );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  lhs_;  //!< Left-hand side dense vector of the addition expression.
   RightOperand rhs_;  //!< Right-hand side sparse vector of the addition expression.
   //**********************************************************************************************

   //**Assignment to dense vectors*****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense vector-sparse vector addition to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side addition expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense vector-sparse
   // vector addition expression to a dense vector.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline void assign( DenseVector<VT,TF>& lhs, const DVecSVecAddExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      if( !IsComputation_v<VT1> && isSame( *lhs, rhs.lhs_ ) ) {
         addAssign( *lhs, rhs.rhs_ );
      }
      else {
         assign   ( *lhs, rhs.lhs_ );
         addAssign( *lhs, rhs.rhs_ );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse vectors****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense vector-sparse vector addition to a sparse vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side addition expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense vector-sparse
   // vector addition expression to a sparse vector.
   */
   template< typename VT >  // Type of the target sparse vector
   friend inline void assign( SparseVector<VT,TF>& lhs, const DVecSVecAddExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      assign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense vectors********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a dense vector-sparse vector addition to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side addition expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense vector-
   // sparse vector addition expression to a dense vector.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline void addAssign( DenseVector<VT,TF>& lhs, const DVecSVecAddExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      addAssign( *lhs, rhs.lhs_ );
      addAssign( *lhs, rhs.rhs_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse vectors*******************************************************
   // No special implementation for the addition assignment to sparse vectors.
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors*****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a dense vector-sparse vector addition to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side addition expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a dense vector-
   // sparse vector addition expression to a dense vector.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline void subAssign( DenseVector<VT,TF>& lhs, const DVecSVecAddExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      subAssign( *lhs, rhs.lhs_ );
      subAssign( *lhs, rhs.rhs_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse vectors****************************************************
   // No special implementation for the subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**Multiplication assignment to dense vectors**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Multiplication assignment of a dense vector-sparse vector addition to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side addition expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a dense
   // vector-sparse vector addition expression to a dense vector.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline void multAssign( DenseVector<VT,TF>& lhs, const DVecSVecAddExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
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

   //**Division assignment to dense vectors********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Division assignment of a dense vector-sparse vector addition to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side addition expression divisor.
   // \return void
   //
   // This function implements the performance optimized division assignment of a dense vector-
   // sparse vector addition expression to a dense vector.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline void divAssign( DenseVector<VT,TF>& lhs, const DVecSVecAddExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      divAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Division assignment to sparse vectors*******************************************************
   // No special implementation for the division assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP assignment to dense vectors*************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a dense vector-sparse vector addition to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side addition expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a dense vector-sparse
   // vector addition expression to a dense vector. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline auto smpAssign( DenseVector<VT,TF>& lhs, const DVecSVecAddExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      if( !IsComputation_v<VT1> && isSame( *lhs, rhs.lhs_ ) ) {
         smpAddAssign( *lhs, rhs.rhs_ );
      }
      else {
         smpAssign   ( *lhs, rhs.lhs_ );
         smpAddAssign( *lhs, rhs.rhs_ );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to sparse vectors************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a dense vector-sparse vector addition to a sparse vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side addition expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a dense vector-sparse
   // vector addition expression to a sparse vector. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename VT >  // Type of the target sparse vector
   friend inline auto smpAssign( SparseVector<VT,TF>& lhs, const DVecSVecAddExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( rhs );
      smpAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to dense vectors****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a dense vector-sparse vector addition to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side addition expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a dense
   // vector-sparse vector addition expression to a dense vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline auto smpAddAssign( DenseVector<VT,TF>& lhs, const DVecSVecAddExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      smpAddAssign( *lhs, rhs.lhs_ );
      smpAddAssign( *lhs, rhs.rhs_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to sparse vectors***************************************************
   // No special implementation for the SMP addition assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP subtraction assignment to dense vectors*************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a dense vector-sparse vector addition to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side addition expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a dense vector-
   // sparse vector addition expression to a dense vector. Due to the explicit application of the
   // SFINAE principle, this function can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline auto smpSubAssign( DenseVector<VT,TF>& lhs, const DVecSVecAddExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      smpSubAssign( *lhs, rhs.lhs_ );
      smpSubAssign( *lhs, rhs.rhs_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to sparse vectors************************************************
   // No special implementation for the SMP subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP multiplication assignment to dense vectors**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP multiplication assignment of a dense vector-sparse vector addition to a dense
   //        vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side addition expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized SMP multiplication assignment of a dense
   // vector-sparse vector addition expression to a dense vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline auto smpMultAssign( DenseVector<VT,TF>& lhs, const DVecSVecAddExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
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

   //**SMP division assignment to dense vectors****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP division assignment of a dense vector-sparse vector addition to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side addition expression divisor.
   // \return void
   //
   // This function implements the performance optimized SMP divison assignment of a dense vector-
   // sparse vector addition expression to a dense vector. Due to the explicit application of the
   // SFINAE principle, this function can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline auto smpDivAssign( DenseVector<VT,TF>& lhs, const DVecSVecAddExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( rhs );
      smpDivAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP division assignment to sparse vectors***************************************************
   // No special implementation for the SMP division assignment to sparse vectors.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( VT1 );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( VT2 );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT1, TF );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT2, TF );
   BLAZE_CONSTRAINT_MUST_FORM_VALID_VECVECADDEXPR( VT1, VT2 );
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
/*!\brief Backend implementation of the addition between a dense vector and a sparse vector
//        (\f$ \vec{a}=\vec{b}+\vec{c} \f$).
// \ingroup dense_vector
//
// \param lhs The left-hand side dense vector for the vector addition.
// \param rhs The right-hand side sparse vector for the vector addition.
// \return The sum of the two vectors.
//
// This function implements a performance optimized treatment of the addition between a dense
// vector and a sparse vector.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side sparse vector
        , bool TF       // Transpose flag
        , DisableIf_t< IsZero_v<VT2> &&
                       IsSame_v< ElementType_t<VT1>, ElementType_t<VT2> > >* = nullptr >
inline const DVecSVecAddExpr<VT1,VT2,TF>
   dvecsvecadd( const DenseVector<VT1,TF>& lhs, const SparseVector<VT2,TF>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).size() == (*rhs).size(), "Invalid vector sizes" );

   return DVecSVecAddExpr<VT1,VT2,TF>( *lhs, *rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the addition between a dense vector and a zero vector
//        (\f$ \vec{a}=\vec{b}+\vec{c} \f$).
// \ingroup dense_vector
//
// \param lhs The left-hand side dense vector for the vector addition.
// \param rhs The right-hand side zero vector for the vector addition.
// \return The sum of the two vectors.
//
// This function implements a performance optimized treatment of the addition between a dense
// vector and a zero vector.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side sparse vector
        , bool TF       // Transpose flag
        , EnableIf_t< IsZero_v<VT2> &&
                      IsSame_v< ElementType_t<VT1>, ElementType_t<VT2> > >* = nullptr >
inline const VT1&
   dvecsvecadd( const DenseVector<VT1,TF>& lhs, const SparseVector<VT2,TF>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( rhs );

   BLAZE_INTERNAL_ASSERT( (*lhs).size() == (*rhs).size(), "Invalid vector sizes" );

   return (*lhs);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition operator for the addition of a dense vector and a sparse vector
//        (\f$ \vec{a}=\vec{b}+\vec{c} \f$).
// \ingroup dense_vector
//
// \param lhs The left-hand side dense vector for the vector addition.
// \param rhs The right-hand side sparse vector for the vector addition.
// \return The sum of the two vectors.
// \exception std::invalid_argument Vector sizes do not match.
//
// This operator represents the addition of a dense vector and a sparse vector:

   \code
   blaze::DynamicVector<double> a, c;
   blaze::CompressedVector<double> b;
   // ... Resizing and initialization
   c = a + b;
   \endcode

// The operator returns an expression representing a dense vector of the higher-order element
// type of the two involved vector element types \a VT1::ElementType and \a VT2::ElementType.
// Both vector types \a VT1 and \a VT2 as well as the two element types \a VT1::ElementType
// and \a VT2::ElementType have to be supported by the AddTrait class template.\n
// In case the current sizes of the two given vectors don't match, a \a std::invalid_argument
// is thrown.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side sparse vector
        , bool TF >     // Transpose flag
inline decltype(auto)
   operator+( const DenseVector<VT1,TF>& lhs, const SparseVector<VT2,TF>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   if( (*lhs).size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   return dvecsvecadd( *lhs, *rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition operator for the addition of a sparse vector and a dense vector
//        (\f$ \vec{a}=\vec{b}+\vec{c} \f$).
// \ingroup dense_vector
//
// \param lhs The left-hand side sparse vector for the vector addition.
// \param rhs The right-hand side dense vector for the vector addition.
// \return The sum of the two vectors.
// \exception std::invalid_argument Vector sizes do not match.
//
// This operator represents the addition of a sparse vector and a dense vector:

   \code
   blaze::CompressedVector<double> a;
   blaze::DynamicVector<double> b, c;
   // ... Resizing and initialization
   c = a + b;
   \endcode

// The operator returns an expression representing a dense vector of the higher-order element
// type of the two involved vector element types \a VT1::ElementType and \a VT2::ElementType.
// Both vector types \a VT1 and \a VT2 as well as the two element types \a VT1::ElementType and
// \a VT2::ElementType have to be supported by the AddTrait class template.\n
// In case the current sizes of the two given vectors don't match, a \a std::invalid_argument
// is thrown.
*/
template< typename VT1  // Type of the left-hand side sparse vector
        , typename VT2  // Type of the right-hand side dense vector
        , bool TF >     // Transpose flag
inline decltype(auto)
   operator+( const SparseVector<VT1,TF>& lhs, const DenseVector<VT2,TF>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   if( (*lhs).size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   return dvecsvecadd( *rhs, *lhs );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition operator for the addition of a dense vector-sparse vector addition
//        expression and a dense vector (\f$ \vec{a}=(\vec{b}+\vec{c})+\vec{d} \f$).
// \ingroup dense_vector
//
// \param lhs The left-hand side dense vector-sparse vector addition.
// \param rhs The right-hand side dense vector.
// \return The sum of the two vectors.
//
// This operator implements a performance optimized treatment of the addition of a dense
// vector-sparse vector addition expression to a dense vector.
*/
template< typename VT1    // Type of the dense vector of the left-hand side expression
        , typename VT2    // Type of the sparse vector of the left-hand side expression
        , bool TF         // Transpose flag of the left-hand side expression
        , typename VT3 >  // Type of the right-hand side dense vector
inline decltype(auto)
   operator+( const DVecSVecAddExpr<VT1,VT2,TF>& lhs, const DenseVector<VT3,TF>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return ( lhs.leftOperand() + (*rhs) ) + lhs.rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction operator for the subtraction of a dense vector-sparse vector addition
//        expression and a dense vector (\f$ \vec{a}=(\vec{b}+\vec{c})-\vec{d} \f$).
// \ingroup dense_vector
//
// \param lhs The left-hand side dense vector-sparse vector addition.
// \param rhs The right-hand side dense vector.
// \return The difference of the two vectors.
//
// This operator implements a performance optimized treatment of the subtraction of a
// dense vector-sparse vector addition expression and a dense vector.
*/
template< typename VT1    // Type of the dense vector of the left-hand side expression
        , typename VT2    // Type of the sparse vector of the left-hand side expression
        , bool TF         // Transpose flag of the left-hand side expression
        , typename VT3 >  // Type of the right-hand side dense vector
inline decltype(auto)
   operator-( const DVecSVecAddExpr<VT1,VT2,TF>& lhs, const DenseVector<VT3,TF>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return ( lhs.leftOperand() - (*rhs) ) + lhs.rightOperand();
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
