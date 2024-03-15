//=================================================================================================
/*!
//  \file blaze/math/expressions/DVecScalarMultExpr.h
//  \brief Header file for the dense vector/scalar multiplication expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DVECSCALARMULTEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DVECSCALARMULTEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <utility>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/Scalar.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/expressions/VecScalarMultExpr.h>
#include <blaze/math/shims/Invert.h>
#include <blaze/math/shims/IsZero.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/SIMD.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/typetraits/HasSIMDMult.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsPadded.h>
#include <blaze/math/typetraits/IsScalar.h>
#include <blaze/math/typetraits/IsTemporary.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/math/typetraits/UnderlyingBuiltin.h>
#include <blaze/math/typetraits/UnderlyingElement.h>
#include <blaze/system/HostDevice.h>
#include <blaze/system/Inline.h>
#include <blaze/system/MacroDisable.h>
#include <blaze/system/Thresholds.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/SameType.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsNumeric.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DVECSCALARMULTEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for dense vector-scalar multiplications.
// \ingroup dense_vector_expression
//
// The DVecScalarMultExpr class represents the compile time expression for multiplications between
// a dense vector and a scalar value.
*/
template< typename VT  // Type of the left-hand side dense vector
        , typename ST  // Type of the right-hand side scalar value
        , bool TF >    // Transpose flag
class DVecScalarMultExpr
   : public VecScalarMultExpr< DenseVector< DVecScalarMultExpr<VT,ST,TF>, TF > >
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   using RT = ResultType_t<VT>;     //!< Result type of the dense vector expression.
   using RN = ReturnType_t<VT>;     //!< Return type of the dense vector expression.
   using ET = ElementType_t<VT>;    //!< Element type of the dense vector expression.
   using CT = CompositeType_t<VT>;  //!< Composite type of the dense vector expression.
   //**********************************************************************************************

   //**Return type evaluation**********************************************************************
   //! Compilation switch for the selection of the subscript operator return type.
   /*! The \a returnExpr compile time constant expression is a compilation switch for the
       selection of the \a ReturnType. If the vector operand returns a temporary vector
       or matrix, \a returnExpr will be set to \a false and the subscript operator will
       return it's result by value. Otherwise \a returnExpr will be set to \a true and
       the subscript operator may return it's result as an expression. */
   static constexpr bool returnExpr = !IsTemporary_v<RN>;

   //! Expression return type for the subscript operator.
   using ExprReturnType = decltype( std::declval<RN>() * std::declval<ST>() );
   //**********************************************************************************************

   //**Serial evaluation strategy******************************************************************
   //! Compilation switch for the serial evaluation strategy of the multiplication expression.
   /*! The \a useAssign compile time constant expression represents a compilation switch for
       the serial evaluation strategy of the multiplication expression. In case the given dense
       vector expression of type \a VT is a computation expression and requires an intermediate
       evaluation, \a useAssign will be set to 1 and the multiplication expression will be
       evaluated via the \a assign function family. Otherwise \a useAssign will be set to 0
       and the expression will be evaluated via the subscript operator. */
   static constexpr bool useAssign = ( IsComputation_v<VT> && RequiresEvaluation_v<VT> );

   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename VT2 >
   static constexpr bool UseAssign_v = useAssign;
   /*! \endcond */
   //**********************************************************************************************

   //**Parallel evaluation strategy****************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! This variable template is a helper for the selection of the parallel evaluation strategy.
       In case either the target vector ir the dense vector operand is not SMP assignable and the
       vector operand is a computation expression that requires an intermediate evaluation, the
       variable is set to 1 and the expression specific evaluation strategy is selected. Otherwise
       the variable is set to 0 and the default strategy is chosen. */
   template< typename VT2 >
   static constexpr bool UseSMPAssign_v =
      ( ( !VT2::smpAssignable || !VT::smpAssignable ) && useAssign );
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this DVecScalarMultExpr instance.
   using This = DVecScalarMultExpr<VT,ST,TF>;

   //! Base type of this DVecScalarMultExpr instance.
   using BaseType = VecScalarMultExpr< DenseVector<This,TF> >;

   using ResultType    = MultTrait_t<RT,ST>;           //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<ResultType>;    //!< Resulting element type.

   //! Return type for expression template evaluations.
   using ReturnType = const If_t< returnExpr, ExprReturnType, ElementType >;

   //! Data type for composite expression templates.
   using CompositeType = If_t< useAssign, const ResultType, const DVecScalarMultExpr& >;

   //! Composite type of the left-hand side dense vector expression.
   using LeftOperand = If_t< IsExpression_v<VT>, const VT, const VT& >;

   //! Composite type of the right-hand side scalar value.
   using RightOperand = ST;
   //**********************************************************************************************

   //**ConstIterator class definition**************************************************************
   /*!\brief Iterator over the elements of the dense vector.
   */
   class ConstIterator
   {
    public:
      //**Type definitions*************************************************************************
      using IteratorCategory = std::random_access_iterator_tag;  //!< The iterator category.
      using ValueType        = ElementType;                      //!< Type of the underlying elements.
      using PointerType      = ElementType*;                     //!< Pointer return type.
      using ReferenceType    = ElementType&;                     //!< Reference return type.
      using DifferenceType   = ptrdiff_t;                        //!< Difference between two iterators.

      // STL iterator requirements
      using iterator_category = IteratorCategory;  //!< The iterator category.
      using value_type        = ValueType;         //!< Type of the underlying elements.
      using pointer           = PointerType;       //!< Pointer return type.
      using reference         = ReferenceType;     //!< Reference return type.
      using difference_type   = DifferenceType;    //!< Difference between two iterators.

      //! ConstIterator type of the dense vector expression.
      using IteratorType = ConstIterator_t<VT>;
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor for the ConstIterator class.
      //
      // \param iterator Iterator to the initial element.
      // \param scalar Scalar of the multiplication expression.
      */
      inline BLAZE_DEVICE_CALLABLE ConstIterator( IteratorType iterator, RightOperand scalar )
         : iterator_( iterator )  // Iterator to the current element
         , scalar_  ( scalar   )  // Scalar of the multiplication expression
      {}
      //*******************************************************************************************

      //**Addition assignment operator*************************************************************
      /*!\brief Addition assignment operator.
      //
      // \param inc The increment of the iterator.
      // \return The incremented iterator.
      */
      inline BLAZE_DEVICE_CALLABLE ConstIterator& operator+=( size_t inc ) {
         iterator_ += inc;
         return *this;
      }
      //*******************************************************************************************

      //**Subtraction assignment operator**********************************************************
      /*!\brief Subtraction assignment operator.
      //
      // \param dec The decrement of the iterator.
      // \return The decremented iterator.
      */
      inline BLAZE_DEVICE_CALLABLE ConstIterator& operator-=( size_t dec ) {
         iterator_ -= dec;
         return *this;
      }
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented iterator.
      */
      inline BLAZE_DEVICE_CALLABLE ConstIterator& operator++() {
         ++iterator_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix increment operator***************************************************************
      /*!\brief Post-increment operator.
      //
      // \return The previous position of the iterator.
      */
      inline BLAZE_DEVICE_CALLABLE const ConstIterator operator++( int ) {
         return ConstIterator( iterator_++, scalar_ );
      }
      //*******************************************************************************************

      //**Prefix decrement operator****************************************************************
      /*!\brief Pre-decrement operator.
      //
      // \return Reference to the decremented iterator.
      */
      inline BLAZE_DEVICE_CALLABLE ConstIterator& operator--() {
         --iterator_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix decrement operator***************************************************************
      /*!\brief Post-decrement operator.
      //
      // \return The previous position of the iterator.
      */
      inline BLAZE_DEVICE_CALLABLE const ConstIterator operator--( int ) {
         return ConstIterator( iterator_--, scalar_ );
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the element at the current iterator position.
      //
      // \return The resulting value.
      */
      inline BLAZE_DEVICE_CALLABLE ReturnType operator*() const {
         return *iterator_ * scalar_;
      }
      //*******************************************************************************************

      //**Load function****************************************************************************
      /*!\brief Access to the SIMD elements of the vector.
      //
      // \return The resulting SIMD element.
      */
      inline auto load() const noexcept {
         return iterator_.load() * set( scalar_ );
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator==( const ConstIterator& rhs ) const {
         return iterator_ == rhs.iterator_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator!=( const ConstIterator& rhs ) const {
         return iterator_ != rhs.iterator_;
      }
      //*******************************************************************************************

      //**Less-than operator***********************************************************************
      /*!\brief Less-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller, \a false if not.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator<( const ConstIterator& rhs ) const {
         return iterator_ < rhs.iterator_;
      }
      //*******************************************************************************************

      //**Greater-than operator********************************************************************
      /*!\brief Greater-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater, \a false if not.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator>( const ConstIterator& rhs ) const {
         return iterator_ > rhs.iterator_;
      }
      //*******************************************************************************************

      //**Less-or-equal-than operator**************************************************************
      /*!\brief Less-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller or equal, \a false if not.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator<=( const ConstIterator& rhs ) const {
         return iterator_ <= rhs.iterator_;
      }
      //*******************************************************************************************

      //**Greater-or-equal-than operator***********************************************************
      /*!\brief Greater-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater or equal, \a false if not.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator>=( const ConstIterator& rhs ) const {
         return iterator_ >= rhs.iterator_;
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two iterators.
      //
      // \param rhs The right-hand side iterator.
      // \return The number of elements between the two iterators.
      */
      inline BLAZE_DEVICE_CALLABLE DifferenceType operator-( const ConstIterator& rhs ) const {
         return iterator_ - rhs.iterator_;
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between a ConstIterator and an integral value.
      //
      // \param it The iterator to be incremented.
      // \param inc The number of elements the iterator is incremented.
      // \return The incremented iterator.
      */
      friend inline BLAZE_DEVICE_CALLABLE const ConstIterator operator+( const ConstIterator& it, size_t inc ) {
         return ConstIterator( it.iterator_ + inc, it.scalar_ );
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between an integral value and a ConstIterator.
      //
      // \param inc The number of elements the iterator is incremented.
      // \param it The iterator to be incremented.
      // \return The incremented iterator.
      */
      friend inline BLAZE_DEVICE_CALLABLE const ConstIterator operator+( size_t inc, const ConstIterator& it ) {
         return ConstIterator( it.iterator_ + inc, it.scalar_ );
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Subtraction between a ConstIterator and an integral value.
      //
      // \param it The iterator to be decremented.
      // \param dec The number of elements the iterator is decremented.
      // \return The decremented iterator.
      */
      friend inline BLAZE_DEVICE_CALLABLE const ConstIterator operator-( const ConstIterator& it, size_t dec ) {
         return ConstIterator( it.iterator_ - dec, it.scalar_ );
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      IteratorType iterator_;  //!< Iterator to the current element.
      RightOperand scalar_;    //!< Scalar of the multiplication expression.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled =
      ( VT::simdEnabled && IsNumeric_v<ET> &&
        ( HasSIMDMult_v<ET,ST> || HasSIMDMult_v<UnderlyingElement_t<ET>,ST> ) );

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = VT::smpAssignable;
   //**********************************************************************************************

   //**SIMD properties*****************************************************************************
   //! The number of elements packed within a single SIMD element.
   static constexpr size_t SIMDSIZE = SIMDTrait<ElementType>::size;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DVecScalarMultExpr class.
   //
   // \param vector The left-hand side dense vector of the multiplication expression.
   // \param scalar The right-hand side scalar of the multiplication expression.
   */
   inline DVecScalarMultExpr( const VT& vector, ST scalar ) noexcept
      : vector_( vector )  // Left-hand side dense vector of the multiplication expression
      , scalar_( scalar )  // Right-hand side scalar of the multiplication expression
   {}
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < vector_.size(), "Invalid vector access index" );
      return vector_[index] * scalar_;
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
      if( index >= vector_.size() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid vector access index" );
      }
      return (*this)[index];
   }
   //**********************************************************************************************

   //**Load function*******************************************************************************
   /*!\brief Access to the SIMD elements of the vector.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return Reference to the accessed values.
   */
   BLAZE_ALWAYS_INLINE auto load( size_t index ) const noexcept {
      BLAZE_INTERNAL_ASSERT( index < vector_.size() , "Invalid vector access index" );
      BLAZE_INTERNAL_ASSERT( index % SIMDSIZE == 0UL, "Invalid vector access index" );
      return vector_.load( index ) * set( scalar_ );
   }
   //**********************************************************************************************

   //**Begin function******************************************************************************
   /*!\brief Returns an iterator to the first non-zero element of the dense vector.
   //
   // \return Iterator to the first non-zero element of the dense vector.
   */
   inline ConstIterator begin() const {
      return ConstIterator( vector_.begin(), scalar_ );
   }
   //**********************************************************************************************

   //**End function********************************************************************************
   /*!\brief Returns an iterator just past the last non-zero element of the dense vector.
   //
   // \return Iterator just past the last non-zero element of the dense vector.
   */
   inline ConstIterator end() const {
      return ConstIterator( vector_.end(), scalar_ );
   }
   //**********************************************************************************************

   //**Size function*******************************************************************************
   /*!\brief Returns the current size/dimension of the vector.
   //
   // \return The size of the vector.
   */
   inline size_t size() const noexcept {
      return vector_.size();
   }
   //**********************************************************************************************

   //**Left operand access*************************************************************************
   /*!\brief Returns the left-hand side dense vector operand.
   //
   // \return The left-hand side dense vector operand.
   */
   inline LeftOperand leftOperand() const noexcept {
      return vector_;
   }
   //**********************************************************************************************

   //**Right operand access************************************************************************
   /*!\brief Returns the right-hand side scalar operand.
   //
   // \return The right-hand side scalar operand.
   */
   inline RightOperand rightOperand() const noexcept {
      return scalar_;
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
      return IsExpression_v<VT> && vector_.canAlias( alias );
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
      return vector_.isAliased( alias );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the operands of the expression are properly aligned in memory.
   //
   // \return \a true in case the operands are aligned, \a false if not.
   */
   inline bool isAligned() const noexcept {
      return vector_.isAligned();
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const noexcept {
      return vector_.canSMPAssign() || ( size() > SMP_DVECSCALARMULT_THRESHOLD );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  vector_;  //!< Left-hand side dense vector of the multiplication expression.
   RightOperand scalar_;  //!< Right-hand side scalar of the multiplication expression.
   //**********************************************************************************************

   //**Assignment to dense vectors*****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense vector-scalar multiplication to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense vector-scalar
   // multiplication expression to a dense vector. Due to the explicit application of the
   // SFINAE principle, this function can only be selected by the compiler in case the vector
   // operand is a computation expression and requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto assign( DenseVector<VT2,TF>& lhs, const DVecScalarMultExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      assign( *lhs, rhs.vector_ );
      assign( *lhs, (*lhs) * rhs.scalar_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse vectors****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense vector-scalar multiplication to a sparse vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense vector-scalar
   // multiplication expression to a sparse vector. Due to the explicit application of the
   // SFINAE principle, this function can only be selected by the compiler in case the vector
   // operand is a computation expression and requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target sparse vector
   friend inline auto assign( SparseVector<VT2,TF>& lhs, const DVecScalarMultExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      assign( *lhs, rhs.vector_ );
      (*lhs) *= rhs.scalar_;
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense vectors********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a dense vector-scalar multiplication to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense vector-
   // scalar multiplication expression to a dense vector. Due to the explicit application of
   // the SFINAE principle, this function can only be selected by the compiler in case the
   // vector operand is a computation expression and requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto addAssign( DenseVector<VT2,TF>& lhs, const DVecScalarMultExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      addAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse vectors*******************************************************
   // No special implementation for the addition assignment to sparse vectors.
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors*****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a dense vector-scalar multiplication to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a dense vector-
   // scalar multiplication expression to a dense vector. Due to the explicit application of the
   // SFINAE principle, this function can only be selected by the compiler in case the vector
   // operand is a computation expression and requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto subAssign( DenseVector<VT2,TF>& lhs, const DVecScalarMultExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      subAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse vectors****************************************************
   // No special implementation for the subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**Multiplication assignment to dense vectors**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Multiplication assignment of a dense vector-scalar multiplication to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a dense
   // vector-scalar multiplication expression to a dense vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case the
   // vector operand is a computation expression and requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto multAssign( DenseVector<VT2,TF>& lhs, const DVecScalarMultExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT2> >
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
   /*!\brief Division assignment of a dense vector-scalar multiplication to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression divisor.
   // \return void
   //
   // This function implements the performance optimized division assignment of a dense vector-
   // scalar multiplication expression to a dense vector. Due to the explicit application of
   // the SFINAE principle, this function can only be selected by the compiler in case the
   // vector operand is a computation expression and requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto divAssign( DenseVector<VT2,TF>& lhs, const DVecScalarMultExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT2> >
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
   /*!\brief SMP assignment of a dense vector-scalar multiplication to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a dense vector-scalar
   // multiplication expression to a dense vector. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto smpAssign( DenseVector<VT2,TF>& lhs, const DVecScalarMultExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      smpAssign( *lhs, rhs.vector_ );
      smpAssign( *lhs, (*lhs) * rhs.scalar_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to sparse vectors************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a dense vector-scalar multiplication to a sparse vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a dense vector-scalar
   // multiplication expression to a sparse vector. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename VT2 >  // Type of the target sparse vector
   friend inline auto smpAssign( SparseVector<VT2,TF>& lhs, const DVecScalarMultExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      smpAssign( *lhs, rhs.vector_ );
      (*lhs) *= rhs.scalar_;
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to dense vectors****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a dense vector-scalar multiplication to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a dense
   // vector-scalar multiplication expression to a dense vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto smpAddAssign( DenseVector<VT2,TF>& lhs, const DVecScalarMultExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( rhs );
      smpAddAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to sparse vectors***************************************************
   // No special implementation for the SMP addition assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP subtraction assignment to dense vectors*************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a dense vector-scalar multiplication to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a dense
   // vector-scalar multiplication expression to a dense vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto smpSubAssign( DenseVector<VT2,TF>& lhs, const DVecScalarMultExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( rhs );
      smpSubAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to sparse vectors************************************************
   // No special implementation for the SMP subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP multiplication assignment to dense vectors**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP multiplication assignment of a dense vector-scalar multiplication to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized SMP multiplication assignment of a dense
   // vector-scalar multiplication expression to a dense vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto smpMultAssign( DenseVector<VT2,TF>& lhs, const DVecScalarMultExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT2> >
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
   /*!\brief SMP division assignment of a dense vector-scalar multiplication to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression divisor.
   // \return void
   //
   // This function implements the performance optimized SMP division assignment of a dense vector-
   // scalar multiplication expression to a dense vector. Due to the explicit application of the
   // SFINAE principle, this function can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto smpDivAssign( DenseVector<VT2,TF>& lhs, const DVecScalarMultExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT2> >
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
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT, TF );
   BLAZE_CONSTRAINT_MUST_BE_SCALAR_TYPE( ST );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ST, RightOperand );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL UNARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Unary minus operator for the negation of a dense vector (\f$ \vec{a} = -\vec{b} \f$).
// \ingroup dense_vector
//
// \param dv The dense vector to be negated.
// \return The negation of the vector.
//
// This operator represents the negation of a dense vector:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = -a;
   \endcode

// The operator returns an expression representing the negation of the given dense vector.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) operator-( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   using ScalarType = UnderlyingBuiltin_t<VT>;
   using ReturnType = const DVecScalarMultExpr<VT,ScalarType,TF>;
   return ReturnType( *dv, ScalarType(-1) );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Multiplication operator for the multiplication of a dense vector and a scalar value
//        (\f$ \vec{a}=\vec{b}*s \f$).
// \ingroup dense_vector
//
// \param vec The left-hand side dense vector for the multiplication.
// \param scalar The right-hand side scalar value for the multiplication.
// \return The scaled result vector.
//
// This operator represents the multiplication between a dense vector and a scalar value:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = a * 1.25;
   \endcode

// The operator returns an expression representing a dense vector of the higher-order element type
// of the involved data types \a VT::ElementType and \a ST. Both data types \a VT::ElementType and
// \a ST have to be supported by the MultTrait class template. Note that this operator only works
// for scalar values of built-in data type.
*/
template< typename VT  // Type of the left-hand side dense vector
        , typename ST  // Type of the right-hand side scalar
        , bool TF      // Transpose flag
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
inline decltype(auto) operator*( const DenseVector<VT,TF>& vec, ST scalar )
{
   BLAZE_FUNCTION_TRACE;

   using ScalarType = MultTrait_t< UnderlyingBuiltin_t<VT>, ST >;
   using ReturnType = const DVecScalarMultExpr<VT,ScalarType,TF>;
   return ReturnType( *vec, scalar );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the multiplication of a scalar value and a dense vector
//        (\f$ \vec{a}=s*\vec{b} \f$).
// \ingroup dense_vector
//
// \param scalar The left-hand side scalar value for the multiplication.
// \param vec The right-hand side vector for the multiplication.
// \return The scaled result vector.
//
// This operator represents the multiplication between a a scalar value and dense vector:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = 1.25 * a;
   \endcode

// The operator returns an expression representing a dense vector of the higher-order element
// type of the involved data types \a ST and \a VT::ElementType. Both data types \a ST and
// \a VT::ElementType have to be supported by the MultTrait class template. Note that this
// operator only works for scalar values of built-in data type.
*/
template< typename ST  // Type of the left-hand side scalar
        , typename VT  // Type of the right-hand side dense vector
        , bool TF      // Transpose flag
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
inline decltype(auto) operator*( ST scalar, const DenseVector<VT,TF>& vec )
{
   BLAZE_FUNCTION_TRACE;

   using ScalarType = MultTrait_t< ST, UnderlyingBuiltin_t<VT> >;
   using ReturnType = const DVecScalarMultExpr<VT,ScalarType,TF>;
   return ReturnType( *vec, scalar );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Normalization of the dense vector (\f$|\vec{a}|=1\f$).
//
// \param vec The given dense vector.
// \return The normalized result vector.
//
// This function represents the normalization of a dense vector:

   \code
   blaze::DynamicVector<double> a;
   // ... Resizing and initialization
   a = normalize( a );
   \endcode

// The function returns an expression representing the normalized dense vector. Note that
// this function only works for floating point vectors. The attempt to use this function for
// an integral vector results in a compile time error.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) normalize( const DenseVector<VT,TF>& vec )
{
   BLAZE_CONSTRAINT_MUST_BE_SCALAR_TYPE( ElementType_t<VT> );

   const auto len ( l2Norm( *vec ) );
   auto ilen( !isZero(len) ? inv(len) : len );

   using ReturnType = const DVecScalarMultExpr<VT,decltype(ilen),TF>;
   return ReturnType( *vec, ilen );
}
//*************************************************************************************************




//=================================================================================================
//
//  ISALIGNED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, typename ST, bool TF >
struct IsAligned< DVecScalarMultExpr<VT,ST,TF> >
   : public IsAligned<VT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISPADDED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, typename ST, bool TF >
struct IsPadded< DVecScalarMultExpr<VT,ST,TF> >
   : public IsPadded<VT>
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
