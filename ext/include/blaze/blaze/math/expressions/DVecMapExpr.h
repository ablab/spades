//=================================================================================================
/*!
//  \file blaze/math/expressions/DVecMapExpr.h
//  \brief Header file for the dense vector map expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DVECMAPEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DVECMAPEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <utility>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/VecMapExpr.h>
#include <blaze/math/Functors.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/SIMD.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/MapTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/typetraits/HasLoad.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsPadded.h>
#include <blaze/math/typetraits/IsPaddingEnabled.h>
#include <blaze/math/typetraits/IsScalar.h>
#include <blaze/math/typetraits/IsSIMDEnabled.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/math/typetraits/UnderlyingScalar.h>
#include <blaze/system/HostDevice.h>
#include <blaze/system/Inline.h>
#include <blaze/system/MacroDisable.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsNumeric.h>
#include <blaze/util/typetraits/IsSame.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DVECMAPEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for the dense vector map() function.
// \ingroup dense_vector_expression
//
// The DVecMapExpr class represents the compile time expression for the evaluation of a custom
// operation on each element of a dense vector via the map() function.
*/
template< typename VT  // Type of the dense vector
        , typename OP  // Type of the custom operation
        , bool TF >    // Transpose flag
class DVecMapExpr
   : public VecMapExpr< DenseVector< DVecMapExpr<VT,OP,TF>, TF > >
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   using RT = ResultType_t<VT>;  //!< Result type of the dense vector expression.
   using ET = ElementType_t<VT>;  //!< Element type of the dense vector expression.
   using RN = ReturnType_t<VT>;  //!< Return type of the dense vector expression.
   //**********************************************************************************************

   //**Serial evaluation strategy******************************************************************
   //! Compilation switch for the serial evaluation strategy of the map expression.
   /*! The \a useAssign compile time constant expression represents a compilation switch for
       the serial evaluation strategy of the map expression. In case the given dense vector
       expression of type \a VT is a computation expression and requires an intermediate
       evaluation, \a useAssign will be set to 1 and the map expression will be evaluated
       via the \a assign function family. Otherwise \a useAssign will be set to 0 and the
       expression will be evaluated via the subscript operator. */
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
       In case either the target vector or the dense vector operand is not SMP assignable and
       the vector operand requires an intermediate evaluation, the variable is set to 1 and the
       expression specific evaluation strategy is selected. Otherwise the variable is set to 0
       and the default strategy is chosen. */
   template< typename VT2 >
   static constexpr bool UseSMPAssign_v =
      ( ( !VT2::smpAssignable || !VT::smpAssignable ) && useAssign );
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this DVecMapExpr instance.
   using This = DVecMapExpr<VT,OP,TF>;

   //! Base type of this DVecMapExpr instance.
   using BaseType = VecMapExpr< DenseVector<This,TF> >;

   using ResultType    = MapTrait_t<RT,OP>;            //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<ResultType>;    //!< Resulting element type.

   //! Return type for expression template evaluations.
   using ReturnType = decltype( std::declval<OP>()( std::declval<RN>() ) );

   //! Data type for composite expression templates.
   using CompositeType = If_t< useAssign, const ResultType, const DVecMapExpr& >;

   //! Composite data type of the dense vector expression.
   using Operand = If_t< IsExpression_v<VT>, const VT, const VT& >;

   //! Data type of the custom unary operation.
   using Operation = OP;
   //**********************************************************************************************

   //**ConstIterator class definition**************************************************************
   /*!\brief Iterator over the elements of the dense vector map expression.
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

      //! ConstIterator type of the left-hand side dense vector expression.
      using IteratorType = ConstIterator_t<VT>;
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor for the ConstIterator class.
      //
      // \param it Iterator to the initial vector element.
      // \param op The custom unary operation.
      */
      inline BLAZE_DEVICE_CALLABLE ConstIterator( IteratorType it, OP op )
         : it_( it )             // Iterator to the current vector element
         , op_( std::move(op) )  // The custom unary operation
      {}
      //*******************************************************************************************

      //**Addition assignment operator*************************************************************
      /*!\brief Addition assignment operator.
      //
      // \param inc The increment of the iterator.
      // \return The incremented iterator.
      */
      inline BLAZE_DEVICE_CALLABLE ConstIterator& operator+=( size_t inc ) {
         it_ += inc;
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
         it_ -= dec;
         return *this;
      }
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented iterator.
      */
      inline BLAZE_DEVICE_CALLABLE ConstIterator& operator++() {
         ++it_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix increment operator***************************************************************
      /*!\brief Post-increment operator.
      //
      // \return The previous position of the iterator.
      */
      inline BLAZE_DEVICE_CALLABLE const ConstIterator operator++( int ) {
         return ConstIterator( it_++, op_ );
      }
      //*******************************************************************************************

      //**Prefix decrement operator****************************************************************
      /*!\brief Pre-decrement operator.
      //
      // \return Reference to the decremented iterator.
      */
      inline BLAZE_DEVICE_CALLABLE ConstIterator& operator--() {
         --it_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix decrement operator***************************************************************
      /*!\brief Post-decrement operator.
      //
      // \return The previous position of the iterator.
      */
      inline BLAZE_DEVICE_CALLABLE const ConstIterator operator--( int ) {
         return ConstIterator( it_--, op_ );
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the element at the current iterator position.
      //
      // \return The resulting value.
      */
      inline BLAZE_DEVICE_CALLABLE ReturnType operator*() const {
         return op_( *it_ );
      }
      //*******************************************************************************************

      //**Load function****************************************************************************
      /*!\brief Access to the SIMD elements of the vector.
      //
      // \return The resulting SIMD element.
      */
      inline auto load() const noexcept {
         return op_.load( it_.load() );
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator==( const ConstIterator& rhs ) const {
         return it_ == rhs.it_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator!=( const ConstIterator& rhs ) const {
         return it_ != rhs.it_;
      }
      //*******************************************************************************************

      //**Less-than operator***********************************************************************
      /*!\brief Less-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller, \a false if not.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator<( const ConstIterator& rhs ) const {
         return it_ < rhs.it_;
      }
      //*******************************************************************************************

      //**Greater-than operator********************************************************************
      /*!\brief Greater-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater, \a false if not.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator>( const ConstIterator& rhs ) const {
         return it_ > rhs.it_;
      }
      //*******************************************************************************************

      //**Less-or-equal-than operator**************************************************************
      /*!\brief Less-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller or equal, \a false if not.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator<=( const ConstIterator& rhs ) const {
         return it_ <= rhs.it_;
      }
      //*******************************************************************************************

      //**Greater-or-equal-than operator***********************************************************
      /*!\brief Greater-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater or equal, \a false if not.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator>=( const ConstIterator& rhs ) const {
         return it_ >= rhs.it_;
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two iterators.
      //
      // \param rhs The right-hand side iterator.
      // \return The number of elements between the two iterators.
      */
      inline BLAZE_DEVICE_CALLABLE DifferenceType operator-( const ConstIterator& rhs ) const {
         return it_ - rhs.it_;
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
         return ConstIterator( it.it_ + inc, it.op_ );
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
         return ConstIterator( it.it_ + inc, it.op_ );
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
         return ConstIterator( it.it_ - dec, it.op_ );
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      IteratorType it_;  //!< Iterator to the current vector element.
      OP           op_;  //!< The custom unary operation.
      //*******************************************************************************************
   };
   //**********************************************************************************************

 public:
   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled =
      ( VT::simdEnabled &&
        If_t< HasSIMDEnabled_v<OP>, GetSIMDEnabled<OP,ET>, HasLoad<OP> >::value );

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = VT::smpAssignable;
   //**********************************************************************************************

   //**SIMD properties*****************************************************************************
   //! The number of elements packed within a single SIMD element.
   static constexpr size_t SIMDSIZE = SIMDTrait<ElementType>::size;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DVecMapExpr class.
   //
   // \param dv The dense vector operand of the map expression.
   // \param op The custom unary operation.
   */
   inline DVecMapExpr( const VT& dv, OP op ) noexcept
      : dv_( dv )             // Dense vector of the map expression
      , op_( std::move(op) )  // The custom unary operation
   {}
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < dv_.size(), "Invalid vector access index" );
      return op_( dv_[index] );
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
      if( index >= dv_.size() ) {
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
      BLAZE_INTERNAL_ASSERT( index < dv_.size()     , "Invalid vector access index" );
      BLAZE_INTERNAL_ASSERT( index % SIMDSIZE == 0UL, "Invalid vector access index" );
      return op_.load( dv_.load( index ) );
   }
   //**********************************************************************************************

   //**Begin function******************************************************************************
   /*!\brief Returns an iterator to the first non-zero element of the dense vector.
   //
   // \return Iterator to the first non-zero element of the dense vector.
   */
   inline ConstIterator begin() const {
      return ConstIterator( dv_.begin(), op_ );
   }
   //**********************************************************************************************

   //**End function********************************************************************************
   /*!\brief Returns an iterator just past the last non-zero element of the dense vector.
   //
   // \return Iterator just past the last non-zero element of the dense vector.
   */
   inline ConstIterator end() const {
      return ConstIterator( dv_.end(), op_ );
   }
   //**********************************************************************************************

   //**Size function*******************************************************************************
   /*!\brief Returns the current size/dimension of the vector.
   //
   // \return The size of the vector.
   */
   inline size_t size() const noexcept {
      return dv_.size();
   }
   //**********************************************************************************************

   //**Operand access******************************************************************************
   /*!\brief Returns the dense vector operand.
   //
   // \return The dense vector operand.
   */
   inline Operand operand() const noexcept {
      return dv_;
   }
   //**********************************************************************************************

   //**Operation access****************************************************************************
   /*!\brief Returns a copy of the custom operation.
   //
   // \return A copy of the custom operation.
   */
   inline Operation operation() const {
      return op_;
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
      return IsExpression_v<VT> && dv_.canAlias( alias );
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
      return dv_.isAliased( alias );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the operands of the expression are properly aligned in memory.
   //
   // \return \a true in case the operands are aligned, \a false if not.
   */
   inline bool isAligned() const noexcept {
      return dv_.isAligned();
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const noexcept {
      return dv_.canSMPAssign();
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   Operand   dv_;  //!< Dense vector of the map expression.
   Operation op_;  //!< The custom unary operation.
   //**********************************************************************************************

   //**Assignment to dense vectors*****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense vector map expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side map expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense vector map
   // expression to a dense vector. Due to the explicit application of the SFINAE principle,
   // this function can only be selected by the compiler in case the operand requires an
   // intermediate evaluation and the underlying numeric data type of the operand and the
   // target vector are identical.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto assign( DenseVector<VT2,TF>& lhs, const DVecMapExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT2> &&
                     IsSame_v< UnderlyingScalar_t<VT>, UnderlyingScalar_t<VT2> > >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      assign( *lhs, rhs.dv_ );
      assign( *lhs, map( *lhs, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to dense vectors*****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense vector map expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side map expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense vector map
   // expression to a dense vector. Due to the explicit application of the SFINAE principle,
   // this function can only be selected by the compiler in case the operand requires an
   // intermediate evaluation and the underlying numeric data type of the operand and the
   // target vector differ.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto assign( DenseVector<VT2,TF>& lhs, const DVecMapExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT2> &&
                     !IsSame_v< UnderlyingScalar_t<VT>, UnderlyingScalar_t<VT2> > >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( RT, TF );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const RT tmp( serial( rhs.dv_ ) );
      assign( *lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse vectors****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense vector map expression to a sparse vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side map expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense vector map
   // expression to a sparse vector. Due to the explicit application of the SFINAE principle,
   // this function can only be selected by the compiler in case the operand requires an
   // intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target sparse vector
   friend inline auto assign( SparseVector<VT2,TF>& lhs, const DVecMapExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( RT, TF );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const RT tmp( serial( rhs.dv_ ) );
      assign( *lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense vectors********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a dense vector map expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side map expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense vector
   // map expression to a dense vector. Due to the explicit application of the SFINAE principle,
   // this function can only be selected by the compiler in case the operand requires an
   // intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto addAssign( DenseVector<VT2,TF>& lhs, const DVecMapExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( RT, TF );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const RT tmp( serial( rhs.dv_ ) );
      addAssign( *lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse vectors*******************************************************
   // No special implementation for the addition assignment to sparse vectors.
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors*****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a dense vector map expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side map expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a dense
   // vector map expression to a dense vector. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the operand
   // requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto subAssign( DenseVector<VT2,TF>& lhs, const DVecMapExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( RT, TF );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const RT tmp( serial( rhs.dv_ ) );
      subAssign( *lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse vectors****************************************************
   // No special implementation for the subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**Multiplication assignment to dense vectors**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Multiplication assignment of a dense vector map expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side map expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a dense
   // vector map expression to a dense vector. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the operand requires
   // an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto multAssign( DenseVector<VT2,TF>& lhs, const DVecMapExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( RT, TF );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const RT tmp( serial( rhs.dv_ ) );
      multAssign( *lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Multiplication assignment to sparse vectors*************************************************
   // No special implementation for the multiplication assignment to sparse vectors.
   //**********************************************************************************************

   //**Division assignment to dense vectors********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Division assignment of a dense vector map expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side map expression divisor.
   // \return void
   //
   // This function implements the performance optimized division assignment of a dense vector
   // map expression to a dense vector. Due to the explicit application of the SFINAE principle,
   // this function can only be selected by the compiler in case the operand requires an
   // intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto divAssign( DenseVector<VT2,TF>& lhs, const DVecMapExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( RT, TF );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const RT tmp( serial( rhs.dv_ ) );
      divAssign( *lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Division assignment to sparse vectors*******************************************************
   // No special implementation for the division assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP assignment to dense vectors*************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a dense vector map expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side map expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a dense vector map
   // expression to a dense vector. Due to the explicit application of the SFINAE principle,
   // this function can only be selected by the compiler in case the expression specific parallel
   // evaluation strategy is selected and the underlying numeric data type of the operand and the
   // target vector are identical
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto smpAssign( DenseVector<VT2,TF>& lhs, const DVecMapExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT2> &&
                     IsSame_v< UnderlyingScalar_t<VT>, UnderlyingScalar_t<VT2> > >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      smpAssign( *lhs, rhs.dv_ );
      smpAssign( *lhs, map( *lhs, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to dense vectors*************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a dense vector map expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side map expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a dense vector map
   // expression to a dense vector. Due to the explicit application of the SFINAE principle,
   // this function can only be selected by the compiler in case the expression specific parallel
   // evaluation strategy is selected and the underlying numeric data type of the operand and the
   // target vector differ.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto smpAssign( DenseVector<VT2,TF>& lhs, const DVecMapExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT2> &&
                     !IsSame_v< UnderlyingScalar_t<VT>, UnderlyingScalar_t<VT2> > >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( RT, TF );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const RT tmp( rhs.dv_ );
      smpAssign( *lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to sparse vectors************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a dense vector map expression to a sparse vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side map expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a dense vector map
   // expression to a sparse vector. Due to the explicit application of the SFINAE principle,
   // this function can only be selected by the compiler in case the expression specific parallel
   // evaluation strategy is selected.
   */
   template< typename VT2 >  // Type of the target sparse vector
   friend inline auto smpAssign( SparseVector<VT2,TF>& lhs, const DVecMapExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( RT, TF );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const RT tmp( rhs.dv_ );
      smpAssign( *lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to dense vectors****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a dense vector map expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side map expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a dense
   // vector map expression to a dense vector. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto smpAddAssign( DenseVector<VT2,TF>& lhs, const DVecMapExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( RT, TF );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const RT tmp( rhs.dv_ );
      smpAddAssign( *lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to sparse vectors***************************************************
   // No special implementation for the SMP addition assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP subtraction assignment to dense vectors*************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a dense vector map expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side map expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a dense
   // vector map expression to a dense vector. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto smpSubAssign( DenseVector<VT2,TF>& lhs, const DVecMapExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( RT, TF );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const RT tmp( rhs.dv_ );
      smpSubAssign( *lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to sparse vectors************************************************
   // No special implementation for the SMP subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP multiplication assignment to dense vectors**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP multiplication assignment of a dense vector map expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side map expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized SMP multiplication assignment of a
   // dense vector map expression to a dense vector. Due to the explicit application of the
   // SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto smpMultAssign( DenseVector<VT2,TF>& lhs, const DVecMapExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( RT, TF );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const RT tmp( rhs.dv_ );
      smpMultAssign( *lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP multiplication assignment to sparse vectors*********************************************
   // No special implementation for the SMP multiplication assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP division assignment to dense vectors****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP division assignment of a dense vector map expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side map expression divisor.
   // \return void
   //
   // This function implements the performance optimized SMP division assignment of a dense
   // vector map expression to a dense vector. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto smpDivAssign( DenseVector<VT2,TF>& lhs, const DVecMapExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( RT, TF );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const RT tmp( rhs.dv_ );
      smpDivAssign( *lhs, map( tmp, rhs.op_ ) );
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
/*!\brief Evaluates the given custom operation on each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector.
// \param op The custom operation.
// \return The custom operation applied to each single element of \a dv.
//
// The \a map() function evaluates the given custom operation on each element of the input
// vector \a dv. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a map() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = map( a, []( double a ){ return std::sqrt( a ); } );
   \endcode
*/
template< typename VT    // Type of the dense vector
        , bool TF        // Transpose flag
        , typename OP >  // Type of the custom operation
inline decltype(auto) map( const DenseVector<VT,TF>& dv, OP op )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const DVecMapExpr<VT,OP,TF>;
   return ReturnType( *dv, std::move(op) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Evaluates the given custom operation on each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector.
// \param op The custom operation.
// \return The custom operation applied to each single element of \a dv.
//
// The \a forEach() function evaluates the given custom operation on each element of the input
// vector \a dv. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a forEach() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = forEach( a, []( double a ){ return std::sqrt( a ); } );
   \endcode
*/
template< typename VT    // Type of the dense vector
        , bool TF        // Transpose flag
        , typename OP >  // Type of the custom operation
inline decltype(auto) forEach( const DenseVector<VT,TF>& dv, OP op )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, std::move(op) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the componentwise minimum of a dense vector \a dv and a scalar.
// \ingroup dense_vector
//
// \param dv The left-hand side dense vector operand.
// \param scalar The right-hand side scalar value.
// \return The resulting dense vector.
//
// This operator computes the componentwise minimum of a dense vector \a dv and a uniform vector
// represented by the scalar value \a scalar. The function returns an expression representing this
// operation.\n
// The following example demonstrates the use of the \a min() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = min( a, 0.0 );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF      // Transpose flag
        , typename ST  // Type of the scalar exponent
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
decltype(auto) min( const DenseVector<VT,TF>& dv, ST scalar )
{
   BLAZE_FUNCTION_TRACE;

   using ET = ElementType_t<VT>;
   using ScalarType = If_t< IsNumeric_v<ET> && IsNumeric_v<ST>, MapTrait_t<ET,ST,Min>, ST >;
   return map( *dv, bind2nd( Min(), ScalarType( scalar ) ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the componentwise minimum of a scalar and a dense vector \a dv.
// \ingroup dense_vector
//
// \param scalar The left-hand side scalar value.
// \param dv The right-hand side dense vector operand.
// \return The resulting dense vector.
//
// This operator computes the componentwise minimum of a uniform vector represented by the scalar
// value \a scalar and a dense vector \a dv. The function returns an expression representing this
// operation.\n
// The following example demonstrates the use of the \a min() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = min( 0.0, a );
   \endcode
*/
template< typename ST  // Type of the scalar exponent
        , typename VT  // Type of the dense vector
        , bool TF      // Transpose flag
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
decltype(auto) min( ST scalar, const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   using ET = ElementType_t<VT>;
   using ScalarType = If_t< IsNumeric_v<ST> && IsNumeric_v<ET>, MapTrait_t<ST,ET,Min>, ST >;
   return map( *dv, bind1st( Min(), ScalarType( scalar ) ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the componentwise maximum of a dense vector \a dv and a scalar.
// \ingroup dense_vector
//
// \param dv The left-hand side dense vector operand.
// \param scalar The right-hand side scalar value.
// \return The resulting dense vector.
//
// This operator computes the componentwise maximum of a dense vector \a dv and a uniform vector
// represented by the scalar value \a scalar. The function returns an expression representing this
// operation.\n
// The following example demonstrates the use of the \a max() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = max( a, 0.0 );
   \endcode

// In case the current sizes of the two given vectors don't match, a \a std::invalid_argument
// is thrown.
*/
template< typename VT  // Type of the dense vector
        , bool TF      // Transpose flag
        , typename ST  // Type of the scalar exponent
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
decltype(auto) max( const DenseVector<VT,TF>& dv, ST scalar )
{
   BLAZE_FUNCTION_TRACE;

   using ET = ElementType_t<VT>;
   using ScalarType = If_t< IsNumeric_v<ET> && IsNumeric_v<ST>, MapTrait_t<ET,ST,Max>, ST >;
   return map( *dv, bind2nd( Max(), ScalarType( scalar ) ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the componentwise maximum of a scalar and a dense vector \a dv.
// \ingroup dense_vector
//
// \param scalar The left-hand side scalar value.
// \param dv The right-hand side dense vector operand.
// \return The resulting dense vector.
//
// This operator computes the componentwise maximum of a uniform vector represented by the scalar
// value \a scalar and a dense vector \a dv. The function returns an expression representing this
// operation.\n
// The following example demonstrates the use of the \a max() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = max( 0.0, a );
   \endcode
*/
template< typename ST  // Type of the scalar exponent
        , typename VT  // Type of the dense vector
        , bool TF      // Transpose flag
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
decltype(auto) max( ST scalar, const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   using ET = ElementType_t<VT>;
   using ScalarType = If_t< IsNumeric_v<ST> && IsNumeric_v<ET>, MapTrait_t<ST,ET,Max>, ST >;
   return map( *dv, bind1st( Max(), ScalarType( scalar ) ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Applies the \a abs() function to each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector.
// \return The resulting dense vector.
//
// This function applies the \a abs() function to each element of the input vector \a dv. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a abs() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = abs( a );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) abs( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, Abs() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Applies the \a sign() function to each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector.
// \return The resulting dense vector.
//
// This function applies the \a sign() function to each element of the input vector \a dv. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a sign() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = sign( a );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) sign( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, Sign() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Applies the \a floor() function to each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector.
// \return The resulting dense vector.
//
// This function applies the \a floor() function to each element of the input vector \a dv. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a floor() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = floor( a );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) floor( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, Floor() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Applies the \a ceil() function to each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector.
// \return The resulting dense vector.
//
// This function applies the \a ceil() function to each element of the input vector \a dv. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a ceil() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = ceil( a );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) ceil( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, Ceil() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Applies the \a trunc() function to each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector.
// \return The resulting dense vector.
//
// This function applies the \a trunc() function to each element of the input vector \a dv. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a trunc() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = trunc( a );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) trunc( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, Trunc() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Applies the \a round() function to each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector.
// \return The resulting dense vector.
//
// This function applies the \a round() function to each element of the input vector \a dv. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a round() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = round( a );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) round( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, Round() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a vector containing the complex conjugate of each single element of \a dv.
// \ingroup dense_vector
//
// \param dv The input vector.
// \return The complex conjugate of each single element of \a dv.
//
// The \a conj() function calculates the complex conjugate of each element of the input vector
// \a dv. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a conj() function:

   \code
   blaze::DynamicVector< complex<double> > a, b;
   // ... Resizing and initialization
   b = conj( a );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) conj( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, Conj() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the conjugate transpose vector of \a dv.
// \ingroup dense_vector
//
// \param dv The input vector.
// \return The conjugate transpose of \a dv.
//
// The \a ctrans() function returns an expression representing the conjugate transpose (also
// called adjoint matrix, Hermitian conjugate matrix or transjugate matrix) of the given input
// vector \a dv.\n
// The following example demonstrates the use of the \a ctrans() function:

   \code
   blaze::DynamicVector< complex<double> > a, b;
   // ... Resizing and initialization
   b = ctrans( a );
   \endcode

// Note that the \a ctrans() function has the same effect as manually applying the \a conj() and
// \a trans function in any order:

   \code
   b = trans( conj( a ) );  // Computing the conjugate transpose vector
   b = conj( trans( a ) );  // Computing the conjugate transpose vector
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) ctrans( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return trans( conj( *dv ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a vector containing the real part of each single element of \a dv.
// \ingroup dense_vector
//
// \param dv The input vector.
// \return The real part of each single element of \a dv.
//
// The \a real() function calculates the real part of each element of the input vector \a dv.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a real() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = real( a );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) real( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, Real() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a vector containing the imaginary part of each single element of \a dv.
// \ingroup dense_vector
//
// \param dv The input vector.
// \return The imaginary part of each single element of \a dv.
//
// The \a imag() function calculates the imaginary part of each element of the input vector \a dv.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a imag() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = imag( a );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) imag( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, Imag() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a vector containing the phase angle of each single element of \a dv.
// \ingroup dense_vector
//
// \param dv The input vector.
// \return The phase angle of each single element of \a dv.
//
// The \a arg() function calculates the phase angle of each element of the input vector \a dv.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a arg() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = arg( a );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) arg( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, Arg() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the square root of each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector; all elements must be in the range \f$[0..\infty)\f$.
// \return The square root of each single element of \a dv.
//
// The \a sqrt() function computes the square root of each element of the input vector \a dv.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a sqrt() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = sqrt( a );
   \endcode

// \note All elements are expected to be in the range \f$[0..\infty)\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) sqrt( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, Sqrt() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the inverse square root of each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector; all elements must be in the range \f$(0..\infty)\f$.
// \return The inverse square root of each single element of \a dv.
//
// The \a invsqrt() function computes the inverse square root of each element of the input vector
// \a dv. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a invsqrt() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = invsqrt( a );
   \endcode

// \note All elements are expected to be in the range \f$(0..\infty)\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) invsqrt( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, InvSqrt() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the cubic root of each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector; all elements must be in the range \f$[0..\infty)\f$.
// \return The cubic root of each single element of \a dv.
//
// The \a cbrt() function computes the cubic root of each element of the input vector \a dv. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a cbrt() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = cbrt( a );
   \endcode

// \note All elements are expected to be in the range \f$[0..\infty)\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) cbrt( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, Cbrt() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the inverse cubic root of each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector; all elements must be in the range \f$(0..\infty)\f$.
// \return The inverse cubic root of each single element of \a dv.
//
// The \a invcbrt() function computes the inverse cubic root of each element of the input vector
// \a dv. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a invcbrt() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = invcbrt( a );
   \endcode

// \note All elements are expected to be in the range \f$(0..\infty)\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) invcbrt( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, InvCbrt() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Restricts each single element of the dense vector \a dv to the range \f$[min..max]\f$.
// \ingroup dense_vector
//
// \param dv The input vector.
// \param min The lower delimiter.
// \param max The upper delimiter.
// \return The vector with restricted elements.
//
// The \a clamp() function restricts each element of the input vector \a dv to the range
// \f$[min..max]\f$. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a clamp() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = clamp( a, -1.0, 1.0 );
   \endcode
*/
template< typename VT    // Type of the dense vector
        , bool TF        // Transpose flag
        , typename DT >  // Type of the delimiters
inline decltype(auto) clamp( const DenseVector<VT,TF>& dv, const DT& min, const DT& max )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, bind2nd( bind3rd( Clamp(), max ), min ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the exponential value for each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector.
// \param exp The scalar exponent.
// \return The exponential value of each single element of \a dv.
//
// The \a pow() function computes the exponential value for each element of the input vector
// \a dv. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a pow() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = pow( a, 4.2 );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF      // Transpose flag
        , typename ST  // Type of the scalar exponent
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
inline decltype(auto) pow( const DenseVector<VT,TF>& dv, ST exp )
{
   BLAZE_FUNCTION_TRACE;

   using ET = ElementType_t<VT>;
   using ScalarType = If_t< IsNumeric_v<ET> && IsNumeric_v<ST>, MultTrait_t<ET,ST>, ST >;
   return map( *dv, blaze::bind2nd( Pow(), ScalarType( exp ) ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes \f$ e^x \f$ for each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector.
// \return The resulting dense vector.
//
// The \a exp() function computes \f$ e^x \f$ for each element of the input vector \a dv. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a exp() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = exp( a );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) exp( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, Exp() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes \f$ 2^x \f$ for each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector.
// \return The resulting dense vector.
//
// The \a exp2() function computes \f$ 2^x \f$ for each element of the input vector \a dv. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a exp2() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = exp2( a );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) exp2( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, Exp2() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes \f$ 10^x \f$ for each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector.
// \return The resulting dense vector.
//
// The \a exp10() function computes \f$ 10^x \f$ for each element of the input vector \a dv. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a exp10() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = exp10( a );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) exp10( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, Exp10() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the natural logarithm for each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector; all elements must be in the range \f$[0..\infty)\f$.
// \return The natural logarithm of each single element of \a dv.
//
// The \a log() function computes natural logarithm for each element of the input vector \a dv.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a log() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = log( a );
   \endcode

// \note All elements are expected to be in the range \f$[0..\infty)\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) log( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, Log() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the binary logarithm for each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector; all elements must be in the range \f$[0..\infty)\f$.
// \return The binary logarithm of each single element of \a dv.
//
// The \a log2() function computes binary logarithm for each element of the input vector \a dv.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a log2() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = log2( a );
   \endcode

// \note All elements are expected to be in the range \f$[0..\infty)\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) log2( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, Log2() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the common logarithm for each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector; all elements must be in the range \f$[0..\infty)\f$.
// \return The common logarithm of each single element of \a dv.
//
// The \a log10() function computes common logarithm for each element of the input vector \a dv.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a log10() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = log10( a );
   \endcode

// \note All elements are expected to be in the range \f$[0..\infty)\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) log10( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, Log10() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the natural logarithm of x+1 for each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector; all elements must be in the range \f$[-1..\infty)\f$.
// \return The natural logarithm of x+1 of each single element of \a dv.
//
// The \a log1p() function computes the natural logarithm of x+1 for each element of the input
// vector \a dv. This may be preferred over the natural logarithm for higher precision computing
// the natural logarithm of a quantity very close to 1. The function returns an expression
// representing this operation.\n
// The following example demonstrates the use of the \a log1p() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = log1p( a );
   \endcode

// \note All elements are expected to be in the range \f$[-1..\infty)\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) log1p( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, Log1p() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the natural logarithm of the absolute value of the gamma function for each
//        single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector; all elements must be in the range \f$[0..\infty)\f$.
// \return The natural logarithm of the absolute value of the gamma function of each single element of \a dv.
//
// The \a lgamma() function computes the natural logarithm of the absolute value of the gamma
// function for each element of the input vector \a dv. The function returns an expression
// representing this operation.\n
// The following example demonstrates the use of the \a lgamma() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = lgamma( a );
   \endcode

// \note All elements are expected to be in the range \f$[0..\infty)\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename VT  // Type of the dense vector
        , bool SO >    // Storage order
inline decltype(auto) lgamma( const DenseVector<VT,SO>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, LGamma() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the sine for each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector.
// \return The sine of each single element of \a dv.
//
// The \a sin() function computes the sine for each element of the input vector \a dv. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a sin() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = sin( a );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) sin( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, Sin() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the inverse sine for each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector; all elements must be in the range \f$[-1..1]\f$.
// \return The inverse sine of each single element of \a dv.
//
// The \a asin() function computes the inverse sine for each element of the input vector \a dv.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a asin() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = asin( a );
   \endcode

// \note All elements are expected to be in the range \f$[-1..1]\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) asin( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, Asin() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the hyperbolic sine for each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector.
// \return The hyperbolic sine of each single element of \a dv.
//
// The \a sinh() function computes the hyperbolic sine for each element of the input vector \a dv.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a sinh() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = sinh( a );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) sinh( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, Sinh() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the inverse hyperbolic sine for each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector.
// \return The inverse hyperbolic sine of each single element of \a dv.
//
// The \a asinh() function computes the inverse hyperbolic sine for each element of the input
// vector \a dv. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a asinh() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = asinh( a );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) asinh( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, Asinh() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the cosine for each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector.
// \return The cosine of each single element of \a dv.
//
// The \a cos() function computes the cosine for each element of the input vector \a dv. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a cos() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = cos( a );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) cos( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, Cos() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the inverse cosine for each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector; all elements must be in the range \f$[-1..1]\f$.
// \return The inverse cosine of each single element of \a dv.
//
// The \a acos() function computes the inverse cosine for each element of the input vector \a dv.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a acos() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = acos( a );
   \endcode

// \note All elements are expected to be in the range \f$[-1..1]\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) acos( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, Acos() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the hyperbolic cosine for each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector.
// \return The hyperbolic cosine of each single element of \a dv.
//
// The \a cosh() function computes the hyperbolic cosine for each element of the input vector
// \a dv. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a cosh() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = cosh( a );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) cosh( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, Cosh() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the inverse hyperbolic cosine for each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector; all elements must be in the range \f$[1..\infty)\f$.
// \return The inverse hyperbolic cosine of each single element of \a dv.
//
// The \a acosh() function computes the inverse hyperbolic cosine for each element of the input
// vector \a dv. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a acosh() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = acosh( a );
   \endcode

// \note All elements are expected to be in the range \f$[1..\infty)\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) acosh( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, Acosh() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the tangent for each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector.
// \return The tangent of each single element of \a dv.
//
// The \a tan() function computes the tangent for each element of the input vector \a dv. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a tan() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = tan( a );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) tan( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, Tan() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the inverse tangent for each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector.
// \return The inverse tangent of each single element of \a dv.
//
// The \a atan() function computes the inverse tangent for each element of the input vector \a dv.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a atan() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = atan( a );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) atan( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, Atan() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the hyperbolic tangent for each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector; all elements must be in the range \f$[-1..1]\f$.
// \return The hyperbolic tangent of each single element of \a dv.
//
// The \a tanh() function computes the hyperbolic tangent for each element of the input vector
// \a dv. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a tanh() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = tanh( a );
   \endcode

// \note All elements are expected to be in the range \f$[-1..1]\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) tanh( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, Tanh() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the inverse hyperbolic tangent for each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector; all elements must be in the range \f$[-1..1]\f$.
// \return The inverse hyperbolic tangent of each single element of \a dv.
//
// The \a atanh() function computes the inverse hyperbolic tangent for each element of the input
// vector \a dv. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a atanh() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = atanh( a );
   \endcode

// \note All elements are expected to be in the range \f$[-1..1]\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) atanh( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, Atanh() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the error function for each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector.
// \return The error function of each single element of \a dv.
//
// The \a erf() function computes the error function for each element of the input vector \a dv.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a erf() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = erf( a );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) erf( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, Erf() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the complementary error function for each single element of the dense vector \a dv.
// \ingroup dense_vector
//
// \param dv The input vector.
// \return The complementary error function of each single element of \a dv.
//
// The \a erfc() function computes the complementary error function for each element of the input
// vector \a dv. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a erfc() function:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = erfc( a );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) erfc( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dv, Erfc() );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Absolute value function for absolute value dense vector expressions.
// \ingroup dense_vector
//
// \param dv The absolute value dense vector expression.
// \return The absolute value of each single element of \a dv.
//
// This function implements a performance optimized treatment of the absolute value operation
// on a dense vector absolute value expression.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) abs( const DVecMapExpr<VT,Abs,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return dv;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Applies the \a sign() function for dense vector \a sign() expressions.
// \ingroup dense_vector
//
// \param dv The dense vector \a sign() expression.
// \return The resulting dense vector.
//
// This function implements a performance optimized treatment of the \a sign() operation on a
// dense vector \a sign() expression.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) sign( const DVecMapExpr<VT,Sign,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return dv;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Applies the \a floor() function to a dense vector \a floor() expressions.
// \ingroup dense_vector
//
// \param dv The dense vector \a floor() expression.
// \return The resulting dense vector.
//
// This function implements a performance optimized treatment of the \a floor() operation on
// a dense vector \a floor() expression.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) floor( const DVecMapExpr<VT,Floor,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return dv;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Applies the \a ceil() function to a dense vector \a ceil() expressions.
// \ingroup dense_vector
//
// \param dv The dense vector \a ceil() expression.
// \return The resulting dense vector.
//
// This function implements a performance optimized treatment of the \a ceil() operation on
// a dense vector \a ceil() expression.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) ceil( const DVecMapExpr<VT,Ceil,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return dv;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Applies the \a trunc() function to a dense vector \a trunc() expressions.
// \ingroup dense_vector
//
// \param dv The dense vector \a trunc() expression.
// \return The resulting dense vector.
//
// This function implements a performance optimized treatment of the \a trunc() operation on
// a dense vector \a trunc() expression.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) trunc( const DVecMapExpr<VT,Trunc,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return dv;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Applies the \a round() function to a dense vector \a round() expressions.
// \ingroup dense_vector
//
// \param dv The dense vector \a round() expression.
// \return The resulting dense vector.
//
// This function implements a performance optimized treatment of the \a round() operation on
// a dense vector \a round() expression.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) round( const DVecMapExpr<VT,Round,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return dv;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Complex conjugate function for complex conjugate dense vector expressions.
// \ingroup dense_vector
//
// \param dv The complex conjugate dense vector expression.
// \return The original dense vector.
//
// This function implements a performance optimized treatment of the complex conjugate operation
// on a dense vector complex conjugate expression. It returns an expression representing the
// original dense vector:

   \code
   blaze::DynamicVector< complex<double> > a, b;
   // ... Resizing and initialization
   b = conj( conj( a ) );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) conj( const DVecMapExpr<VT,Conj,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return dv.operand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Complex conjugate function for conjugate transpose dense vector expressions.
// \ingroup dense_vector
//
// \param dv The conjugate transpose dense vector expression.
// \return The transpose dense vector.
//
// This function implements a performance optimized treatment of the complex conjugate operation
// on a dense vector conjugate transpose expression. It returns an expression representing the
// transpose of the dense vector:

   \code
   blaze::DynamicVector< complex<double> > a, b;
   // ... Resizing and initialization
   b = conj( ctrans( a ) );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) conj( const DVecTransExpr<DVecMapExpr<VT,Conj,TF>,!TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return trans( dv.operand().operand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Real part function for real part dense vector expressions.
// \ingroup dense_vector
//
// \param dv The real part dense vector expression.
// \return The real part of each single element of \a dv.
//
// This function implements a performance optimized treatment of the real part operation on
// a dense vector real part expression.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) real( const DVecMapExpr<VT,Real,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return dv;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Addition operator for the addition of a dense vector and a scalar value
//        (\f$ \vec{a}=\vec{b}+s \f$).
// \ingroup dense_vector
//
// \param vec The left-hand side dense vector for the addition.
// \param scalar The right-hand side scalar value for the addition.
// \return The vector sum.
//
// This operator represents the elementwise addition of a dense vector and a uniform vector
// represented by a scalar value:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = a + 1.25;
   \endcode

// The operator returns an expression representing a dense vector of the higher-order element type
// of the involved data types \a VT::ElementType and \a ST. Both data types \a VT::ElementType and
// \a ST have to be supported by the AddTrait class template. Note that this operator only works
// for scalar values of built-in data type.
*/
template< typename VT  // Type of the left-hand side dense vector
        , bool TF      // Transpose flag of the left-hand side dense vector
        , typename ST  // Type of the right-hand side scalar
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
inline decltype(auto) operator+( const DenseVector<VT,TF>& vec, ST scalar )
{
   BLAZE_FUNCTION_TRACE;

   using ET = ElementType_t<VT>;
   using ScalarType = If_t< IsNumeric_v<ET> && IsNumeric_v<ST>, AddTrait_t<ET,ST>, ST >;
   return map( *vec, blaze::bind2nd( Add{}, ScalarType( scalar ) ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition operator for the addition of a scalar value and a dense vector
//        (\f$ \vec{a}=s+\vec{b} \f$).
// \ingroup dense_vector
//
// \param scalar The left-hand side scalar value for the addition.
// \param vec The right-hand side dense vector for the addition.
// \return The vector sum.
//
// This operator represents the elementwise addition of a uniform vector represented by a scalar
// value and a dense vector:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = 1.25 + a;
   \endcode

// The operator returns an expression representing a dense vector of the higher-order element type
// of the involved data types \a VT::ElementType and \a ST. Both data types \a VT::ElementType and
// \a ST have to be supported by the AddTrait class template. Note that this operator only works
// for scalar values of built-in data type.
*/
template< typename ST  // Type of the left-hand side scalar
        , typename VT  // Type of the right-hand side dense vector
        , bool TF      // Transpose flag of the right-hand side dense vector
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
inline decltype(auto) operator+( ST scalar, const DenseVector<VT,TF>& vec )
{
   BLAZE_FUNCTION_TRACE;

   using ET = ElementType_t<VT>;
   using ScalarType = If_t< IsNumeric_v<ST> && IsNumeric_v<ET>, AddTrait_t<ST,ET>, ST >;
   return map( *vec, blaze::bind1st( Add{}, ScalarType( scalar ) ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction operator for the subtraction of a dense vector and a scalar value
//        (\f$ \vec{a}=\vec{b}-s \f$).
// \ingroup dense_vector
//
// \param vec The left-hand side dense vector for the subtraction.
// \param scalar The right-hand side scalar value for the subtraction.
// \return The vector difference.
//
// This operator represents the elementwise subtraction of a uniform vector represented by a
// scalar value from a dense vector:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = a - 1.25;
   \endcode

// The operator returns an expression representing a dense vector of the higher-order element type
// of the involved data types \a VT::ElementType and \a ST. Both data types \a VT::ElementType and
// \a ST have to be supported by the SubTrait class template. Note that this operator only works
// for scalar values of built-in data type.
*/
template< typename VT  // Type of the left-hand side dense vector
        , bool TF      // Transpose flag of the left-hand side dense vector
        , typename ST  // Type of the right-hand side scalar
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
inline decltype(auto) operator-( const DenseVector<VT,TF>& vec, ST scalar )
{
   BLAZE_FUNCTION_TRACE;

   using ET = ElementType_t<VT>;
   using ScalarType = If_t< IsNumeric_v<ET> && IsNumeric_v<ST>, SubTrait_t<ET,ST>, ST >;
   return map( *vec, blaze::bind2nd( Sub{}, ScalarType( scalar ) ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction operator for the subtraction of a scalar value and a dense vector
//        (\f$ \vec{a}=s-\vec{b} \f$).
// \ingroup dense_vector
//
// \param scalar The left-hand side scalar value for the subtraction.
// \param vec The right-hand side dense vector for the subtraction.
// \return The vector difference.
//
// This operator represents the elementwise subtraction of a dense vector from a uniform vector
// represented by a scalar value:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = 1.25 - a;
   \endcode

// The operator returns an expression representing a dense vector of the higher-order element type
// of the involved data types \a VT::ElementType and \a ST. Both data types \a VT::ElementType and
// \a ST have to be supported by the SubTrait class template. Note that this operator only works
// for scalar values of built-in data type.
*/
template< typename ST  // Type of the left-hand side scalar
        , typename VT  // Type of the right-hand side dense vector
        , bool TF      // Transpose flag of the right-hand side dense vector
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
inline decltype(auto) operator-( ST scalar, const DenseVector<VT,TF>& vec )
{
   BLAZE_FUNCTION_TRACE;

   using ET = ElementType_t<VT>;
   using ScalarType = If_t< IsNumeric_v<ST> && IsNumeric_v<ET>, SubTrait_t<ST,ET>, ST >;
   return map( *vec, blaze::bind1st( Sub{}, ScalarType( scalar ) ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division operator for the division of a scalar value and a dense vector
//        (\f$ \vec{a}=s/\vec{b} \f$).
// \ingroup dense_vector
//
// \param scalar The left-hand side scalar value for the division.
// \param vec The right-hand side dense vector for the division.
// \return The vector quotient.
//
// This operator represents the elementwise division of a uniform vector represented by a scalar
// value and a dense vector:

   \code
   blaze::DynamicVector<double> a, b;
   // ... Resizing and initialization
   b = 1.25 / a;
   \endcode

// The operator returns an expression representing a dense vector of the higher-order element type
// of the involved data types \a VT::ElementType and \a ST. Both data types \a VT::ElementType and
// \a ST have to be supported by the DivTrait class template. Note that this operator only works
// for scalar values of built-in data type.
*/
template< typename ST  // Type of the left-hand side scalar
        , typename VT  // Type of the right-hand side dense vector
        , bool TF      // Transpose flag of the right-hand side dense vector
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
inline decltype(auto) operator/( ST scalar, const DenseVector<VT,TF>& vec )
{
   BLAZE_FUNCTION_TRACE;

   using ET = ElementType_t<VT>;
   using ScalarType = If_t< IsNumeric_v<ST> && IsNumeric_v<ET>, DivTrait_t<ST,ET>, ST >;
   return map( *vec, blaze::bind1st( Div{}, ScalarType( scalar ) ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Left-shift operator for the uniform left-shift of a dense vector.
// \ingroup dense_vector
//
// \param vec The dense vector for the uniform left-shift operation.
// \param count The number of bits to shift all vector elements.
// \return The resulting vector.
//
// This operator represents the uniform left-shift of all elements of a dense vector:

   \code
   blaze::DynamicVector<unsigned int> a, b;
   // ... Resizing and initialization
   b = a << 3;
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) operator<<( const DenseVector<VT,TF>& vec, int count )
{
   BLAZE_FUNCTION_TRACE;

   return map( *vec, ShiftLI( count ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Right-shift operator for the uniform right-shift of a dense vector.
// \ingroup dense_vector
//
// \param vec The dense vector for the uniform right-shift operation.
// \param count The number of bits to shift all vector elements.
// \return The resulting vector.
//
// This operator represents the uniform right-shift of all elements of a dense vector:

   \code
   blaze::DynamicVector<unsigned int> a, b;
   // ... Resizing and initialization
   b = a >> 3;
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) operator>>( const DenseVector<VT,TF>& vec, int count )
{
   BLAZE_FUNCTION_TRACE;

   return map( *vec, ShiftRI( count ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise AND operator for the bitwise AND of a dense vector and a scalar value.
// \ingroup dense_vector
//
// \param vec The left-hand side dense vector for the bitwise AND.
// \param scalar The right-hand side scalar value for the bitwise AND.
// \return The resulting vector.
//
// This operator represents the bitwise AND of a scalar value with all elements of a dense vector:

   \code
   blaze::DynamicVector<unsigned int> a, b;
   // ... Resizing and initialization
   b = a & 7U;
   \endcode
*/
template< typename VT  // Type of the left-hand side dense vector
        , bool TF      // Transpose flag
        , typename ST  // Type of the right-hand side scalar
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
inline decltype(auto) operator&( const DenseVector<VT,TF>& vec, ST scalar )
{
   BLAZE_FUNCTION_TRACE;

   return map( *vec, blaze::bind2nd( Bitand{}, scalar ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise OR operator for the bitwise OR of a dense vector and a scalar value.
// \ingroup dense_vector
//
// \param vec The left-hand side dense vector for the bitwise OR.
// \param scalar The right-hand side scalar value for the bitwise OR.
// \return The resulting vector.
//
// This operator represents the bitwise OR of a scalar value with all elements of a dense vector:

   \code
   blaze::DynamicVector<unsigned int> a, b;
   // ... Resizing and initialization
   b = a | 7U;
   \endcode
*/
template< typename VT  // Type of the left-hand side dense vector
        , bool TF      // Transpose flag
        , typename ST  // Type of the right-hand side scalar
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
inline decltype(auto) operator|( const DenseVector<VT,TF>& vec, ST scalar )
{
   BLAZE_FUNCTION_TRACE;

   return map( *vec, blaze::bind2nd( Bitor{}, scalar ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise XOR operator for the bitwise XOR of a dense vector and a scalar value.
// \ingroup dense_vector
//
// \param vec The left-hand side dense vector for the bitwise XOR.
// \param scalar The right-hand side scalar value for the bitwise XOR.
// \return The resulting vector.
//
// This operator represents the bitwise XOR of a scalar value with all elements of a dense vector:

   \code
   blaze::DynamicVector<unsigned int> a, b;
   // ... Resizing and initialization
   b = a ^ 7U;
   \endcode
*/
template< typename VT  // Type of the left-hand side dense vector
        , bool TF      // Transpose flag
        , typename ST  // Type of the right-hand side scalar
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
inline decltype(auto) operator^( const DenseVector<VT,TF>& vec, ST scalar )
{
   BLAZE_FUNCTION_TRACE;

   return map( *vec, blaze::bind2nd( Bitxor{}, scalar ) );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL LOGICAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Logical NOT operator for the logical NOT of a dense vector.
// \ingroup dense_vector
//
// \param vec The dense vector for the logical NOT.
// \return The negated vector.
//
// This operator represents the logical NOT of all elements of a dense vector:

   \code
   blaze::DynamicVector<bool> a, b;
   // ... Resizing and initialization
   b = !a;
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) operator!( const DenseVector<VT,TF>& vec )
{
   BLAZE_FUNCTION_TRACE;

   return map( *vec, Not{} );
}
//*************************************************************************************************




//=================================================================================================
//
//  ISALIGNED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, typename OP, bool TF >
struct IsAligned< DVecMapExpr<VT,OP,TF> >
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
template< typename VT, typename OP, bool TF >
struct IsPadded< DVecMapExpr<VT,OP,TF> >
   : public BoolConstant< IsPadded_v<VT> && IsPaddingEnabled_v<OP> >
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
