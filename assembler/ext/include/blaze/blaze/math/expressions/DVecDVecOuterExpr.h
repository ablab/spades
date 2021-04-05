//=================================================================================================
/*!
//  \file blaze/math/expressions/DVecDVecOuterExpr.h
//  \brief Header file for the dense vector/dense vector outer map expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DVECDVECOUTEREXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DVECDVECOUTEREXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <utility>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/ColumnVector.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
#include <blaze/math/constraints/RowVector.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/constraints/Symmetric.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/VecTVecMapExpr.h>
#include <blaze/math/functors/Add.h>
#include <blaze/math/functors/Div.h>
#include <blaze/math/functors/Mult.h>
#include <blaze/math/functors/Sub.h>
#include <blaze/math/shims/PrevMultiple.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/SIMD.h>
#include <blaze/math/traits/MapTrait.h>
#include <blaze/math/typetraits/HasLoad.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsColumnMajorMatrix.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsPadded.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsSIMDCombinable.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/math/typetraits/Size.h>
#include <blaze/system/HostDevice.h>
#include <blaze/system/Inline.h>
#include <blaze/system/MacroDisable.h>
#include <blaze/system/Optimizations.h>
#include <blaze/system/Thresholds.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsSame.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DVECDVECMAPEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for the dense vector-dense vector outer map() function.
// \ingroup dense_matrix_expression
//
// The DVecDVecOuterExpr class represents the compile time expression for the pairwise (outer)
// evaluation of a binary custom operation on the elements of two dense vectors via the map()
// function.
*/
template< typename VT1   // Type of the left-hand side dense vector
        , typename VT2   // Type of the right-hand side dense vector
        , typename OP >  // Type of the custom operation
class DVecDVecOuterExpr
   : public VecTVecMapExpr< DenseMatrix< DVecDVecOuterExpr<VT1,VT2,OP>, false > >
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   using RT1 = ResultType_t<VT1>;     //!< Result type of the left-hand side dense vector expression.
   using RT2 = ResultType_t<VT2>;     //!< Result type of the right-hand side dense vector expression.
   using ET1 = ElementType_t<VT1>;    //!< Element type of the left-hand side dense vector expression.
   using ET2 = ElementType_t<VT2>;    //!< Element type of the right-hand side dense vector expression.
   using RN1 = ReturnType_t<VT1>;     //!< Return type of the left-hand side dense vector expression.
   using RN2 = ReturnType_t<VT2>;     //!< Return type of the right-hand side dense vector expression.
   using CT1 = CompositeType_t<VT1>;  //!< Composite type of the left-hand side dense vector expression.
   using CT2 = CompositeType_t<VT2>;  //!< Composite type of the right-hand side dense vector expression.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the left-hand side dense vector expression.
   static constexpr bool evaluateLeft = ( IsComputation_v<VT1> || RequiresEvaluation_v<VT1> );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the right-hand side dense vector expression.
   static constexpr bool evaluateRight = ( IsComputation_v<VT2> || RequiresEvaluation_v<VT2> );
   //**********************************************************************************************

   //**Serial evaluation strategy******************************************************************
   //! Compilation switch for the serial evaluation strategy of the map expression.
   /*! The \a useAssign compile time constant expression represents a compilation switch for the
       serial evaluation strategy of the map expression. In case either of the two dense vector
       operands requires an intermediate evaluation, \a useAssign will be set to \a true and the
       map expression will be evaluated via the \a assign function family. Otherwise \a useAssign
       will be set to \a false and the expression will be evaluated via the function call operator. */
   static constexpr bool useAssign = ( evaluateLeft || evaluateRight );
   //**********************************************************************************************

   //**Parallel evaluation strategy****************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! This variable template is a helper for the selection of the parallel evaluation strategy.
       In case the right-hand side vector operand requires an intermediate evaluation, the variable
       will be set to \a true, otherwise it will be \a false. */
   template< typename VT >
   static constexpr bool UseSMPAssign_v = evaluateRight;
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! In case all three involved data types are suited for a vectorized computation of the
       outer product, the variable will be set to \a true, otherwise it will be \a false. */
   template< typename T1, typename T2, typename T3 >
   static constexpr bool UseVectorizedKernel_v =
      ( useOptimizedKernels &&
        T1::simdEnabled && T2::simdEnabled && T3::simdEnabled &&
        IsSIMDCombinable_v< ElementType_t<T1>
                          , ElementType_t<T2>
                          , ElementType_t<T3> > &&
        If_t< HasSIMDEnabled_v<OP>, GetSIMDEnabled<OP,ET1,ET2>, HasLoad<OP> >::value );
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! In case no vectorized computation is possible, the variable will be set to \a true,
       otherwise it will be \a false. */
   template< typename T1, typename T2, typename T3 >
   static constexpr bool UseDefaultKernel_v = !UseVectorizedKernel_v<T1,T2,T3>;
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this DVecDVecOuterExpr instance.
   using This = DVecDVecOuterExpr<VT1,VT2,OP>;

   //! Base type of this DVecDVecOuterExpr instance.
   using BaseType = VecTVecMapExpr< DenseMatrix<This,false> >;

   using ResultType    = MapTrait_t<RT1,RT2,OP>;         //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<ResultType>;     //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;    //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<ResultType>;      //!< Resulting element type.

   //! Return type for expression template evaluations.
   using ReturnType = decltype( std::declval<OP>()( std::declval<RN1>(), std::declval<RN2>() ) );

   //! Data type for composite expression templates.
   using CompositeType = If_t< useAssign, const ResultType, const DVecDVecOuterExpr& >;

   //! Composite type of the left-hand side dense vector expression.
   using LeftOperand = If_t< IsExpression_v<VT1>, const VT1, const VT1& >;

   //! Composite type of the right-hand side dense vector expression.
   using RightOperand = If_t< IsExpression_v<VT2>, const VT2, const VT2& >;

   //! Data type of the custom binary operation.
   using Operation = OP;

   //! Type for the assignment of the left-hand side dense vector operand.
   using LT = If_t< RequiresEvaluation_v<VT1>, const RT1, CT1 >;

   //! Type for the assignment of the right-hand side dense vector operand.
   using RT = If_t< RequiresEvaluation_v<VT2>, const RT2, CT2 >;
   //**********************************************************************************************

   //**ConstIterator class definition**************************************************************
   /*!\brief Iterator over the elements of the dense vector outer map expression.
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
      using LeftIteratorType = ConstIterator_t<VT1>;

      //! ConstIterator type of the right-hand side dense vector expression.
      using RightIteratorType = ConstIterator_t<VT2>;
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor for the ConstIterator class.
      //
      // \param left Iterator to the initial left-hand side element.
      // \param right Iterator to the initial right-hand side element.
      // \param op The custom binary operation.
      */
      inline ConstIterator( LeftIteratorType left, RightIteratorType right, OP op )
         : left_ ( left  )          // Iterator to the current left-hand side element
         , right_( right )          // Iterator to the current right-hand side element
         , op_   ( std::move(op) )  // The custom binary operation
      {}
      //*******************************************************************************************

      //**Addition assignment operator*************************************************************
      /*!\brief Addition assignment operator.
      //
      // \param inc The increment of the iterator.
      // \return The incremented iterator.
      */
      inline BLAZE_DEVICE_CALLABLE ConstIterator& operator+=( size_t inc ) {
         right_ += inc;
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
         right_ -= dec;
         return *this;
      }
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented iterator.
      */
      inline BLAZE_DEVICE_CALLABLE ConstIterator& operator++() {
         ++right_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix increment operator***************************************************************
      /*!\brief Post-increment operator.
      //
      // \return The previous position of the iterator.
      */
      inline BLAZE_DEVICE_CALLABLE const ConstIterator operator++( int ) {
         return ConstIterator( left_, right_++, op_ );
      }
      //*******************************************************************************************

      //**Prefix decrement operator****************************************************************
      /*!\brief Pre-decrement operator.
      //
      // \return Reference to the decremented iterator.
      */
      inline BLAZE_DEVICE_CALLABLE ConstIterator& operator--() {
         --right_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix decrement operator***************************************************************
      /*!\brief Post-decrement operator.
      //
      // \return The previous position of the iterator.
      */
      inline BLAZE_DEVICE_CALLABLE const ConstIterator operator--( int ) {
         return ConstIterator( left_, right_--, op_ );
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the element at the current iterator position.
      //
      // \return The resulting value.
      */
      inline ReturnType operator*() const {
         return op_( *left_, *right_ );
      }
      //*******************************************************************************************

      //**Load function****************************************************************************
      /*!\brief Access to the SIMD elements of the vector.
      //
      // \return The resulting SIMD element.
      */
      inline auto load() const noexcept {
         return op_.load( set( *left_ ), right_.load() );
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      inline bool operator==( const ConstIterator& rhs ) const {
         return right_ == rhs.right_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      inline bool operator!=( const ConstIterator& rhs ) const {
         return right_ != rhs.right_;
      }
      //*******************************************************************************************

      //**Less-than operator***********************************************************************
      /*!\brief Less-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller, \a false if not.
      */
      inline bool operator<( const ConstIterator& rhs ) const {
         return right_ < rhs.right_;
      }
      //*******************************************************************************************

      //**Greater-than operator********************************************************************
      /*!\brief Greater-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater, \a false if not.
      */
      inline bool operator>( const ConstIterator& rhs ) const {
         return right_ > rhs.right_;
      }
      //*******************************************************************************************

      //**Less-or-equal-than operator**************************************************************
      /*!\brief Less-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller or equal, \a false if not.
      */
      inline bool operator<=( const ConstIterator& rhs ) const {
         return right_ <= rhs.right_;
      }
      //*******************************************************************************************

      //**Greater-or-equal-than operator***********************************************************
      /*!\brief Greater-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater or equal, \a false if not.
      */
      inline bool operator>=( const ConstIterator& rhs ) const {
         return right_ >= rhs.right_;
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two iterators.
      //
      // \param rhs The right-hand side iterator.
      // \return The number of elements between the two iterators.
      */
      inline DifferenceType operator-( const ConstIterator& rhs ) const {
         return right_ - rhs.right_;
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between a ConstIterator and an integral value.
      //
      // \param it The iterator to be incremented.
      // \param inc The number of elements the iterator is incremented.
      // \return The incremented iterator.
      */
      friend inline const ConstIterator operator+( const ConstIterator& it, size_t inc ) {
         return ConstIterator( it.left_ + inc, it.right_ + inc, it.op_ );
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between an integral value and a ConstIterator.
      //
      // \param inc The number of elements the iterator is incremented.
      // \param it The iterator to be incremented.
      // \return The incremented iterator.
      */
      friend inline const ConstIterator operator+( size_t inc, const ConstIterator& it ) {
         return ConstIterator( it.left_, it.right_ + inc, it.op_ );
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Subtraction between a ConstIterator and an integral value.
      //
      // \param it The iterator to be decremented.
      // \param dec The number of elements the iterator is decremented.
      // \return The decremented iterator.
      */
      friend inline const ConstIterator operator-( const ConstIterator& it, size_t dec ) {
         return ConstIterator( it.left_, it.right_ - dec, it.op_ );
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      LeftIteratorType  left_;   //!< Iterator to the current left-hand side element.
      RightIteratorType right_;  //!< Iterator to the current right-hand side element.
      OP                op_;     //!< The custom binary operation.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled =
      ( VT1::simdEnabled && VT2::simdEnabled &&
        If_t< HasSIMDEnabled_v<OP>, GetSIMDEnabled<OP,ET1,ET2>, HasLoad<OP> >::value );

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = ( VT1::smpAssignable && VT2::smpAssignable );
   //**********************************************************************************************

   //**SIMD properties*****************************************************************************
   //! The number of elements packed within a single SIMD element.
   static constexpr size_t SIMDSIZE = SIMDTrait<ElementType>::size;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DVecDVecOuterExpr class.
   //
   // \param lhs The left-hand side dense vector operand of the map expression.
   // \param rhs The right-hand side dense vector operand of the map expression.
   // \param op The custom binary operation.
   */
   inline DVecDVecOuterExpr( const VT1& lhs, const VT2& rhs, OP op ) noexcept
      : lhs_( lhs )            // Left-hand side dense vector of the map expression
      , rhs_( rhs )            // Right-hand side dense vector of the map expression
      , op_ ( std::move(op) )  // The custom binary operation
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

      return op_( lhs_[i], rhs_[j] );
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

   //**Load function*******************************************************************************
   /*!\brief Access to the SIMD elements of the matrix.
   //
   // \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
   // \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
   // \return Reference to the accessed values.
   */
   BLAZE_ALWAYS_INLINE auto load( size_t i, size_t j ) const noexcept {
      BLAZE_INTERNAL_ASSERT( i < lhs_.size()    , "Invalid row access index"    );
      BLAZE_INTERNAL_ASSERT( j < rhs_.size()    , "Invalid column access index" );
      BLAZE_INTERNAL_ASSERT( j % SIMDSIZE == 0UL, "Invalid column access index" );
      return op_.load( set( lhs_[i] ), rhs_.load( j ) );
   }
   //**********************************************************************************************

   //**Begin function******************************************************************************
   /*!\brief Returns an iterator to the first non-zero element of row \a i.
   //
   // \param i The row index.
   // \return Iterator to the first non-zero element of row \a i.
   */
   inline ConstIterator begin( size_t i ) const {
      BLAZE_INTERNAL_ASSERT( i < lhs_.size(), "Invalid row access index" );
      return ConstIterator( lhs_.begin()+i, rhs_.begin(), op_ );
   }
   //**********************************************************************************************

   //**End function********************************************************************************
   /*!\brief Returns an iterator just past the last non-zero element of row \a i.
   //
   // \param i The row index.
   // \return Iterator just past the last non-zero element of row \a i.
   */
   inline ConstIterator end( size_t i ) const {
      BLAZE_INTERNAL_ASSERT( i < lhs_.size(), "Invalid row access index" );
      return ConstIterator( lhs_.begin()+i, rhs_.end(), op_ );
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
   /*!\brief Returns the right-hand side dense vector operand.
   //
   // \return The right-hand side dense vector operand.
   */
   inline RightOperand rightOperand() const noexcept {
      return rhs_;
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
      return lhs_.isAligned() && rhs_.isAligned();
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const noexcept {
      return ( rows() * columns() >= SMP_DVECDVECOUTER_THRESHOLD );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  lhs_;  //!< Left-hand side dense vector of the map expression.
   RightOperand rhs_;  //!< Right-hand side dense vector of the map expression.
   Operation    op_;   //!< The custom binary operation.
   //**********************************************************************************************

   //**Assignment to dense matrices****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense vector-dense vector outer map expression to a row-major dense
   //        matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side outer map expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense vector-dense
   // vector outer map expression to a row-major dense matrix. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case either
   // of the two operands requires an intermediate evaluation.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline auto assign( DenseMatrix<MT,false>& lhs, const DVecDVecOuterExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT x( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense vector operand
      RT y( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (*lhs).rows()   , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == (*lhs).columns(), "Invalid vector size" );

      DVecDVecOuterExpr::selectAssignKernel( *lhs, x, y, rhs.op_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default assignment to row-major dense matrices**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default assignment of a dense vector-dense vector outer map expression to a row-major
   //        dense matrix.
   // \ingroup dense_matrix
   //
   // \param A The target left-hand side dense matrix.
   // \param x The left-hand side dense vector operand.
   // \param y The right-hand side dense vector operand.
   // \param op The custom binary operation.
   // \return void
   //
   // This function implements the default assignment kernel for the dense vector-dense vector
   // outer map expression.
   */
   template< typename MT     // Type of the left-hand side target matrix
           , typename VT3    // Type of the left-hand side vector operand
           , typename VT4 >  // Type of the right-hand side vector operand
   static inline auto selectAssignKernel( MT& A, const VT3& x, const VT4& y, OP op )
      -> EnableIf_t< IsRowMajorMatrix_v<MT> && UseDefaultKernel_v<MT,VT3,VT4> >
   {
      const size_t M( A.rows() );
      const size_t N( A.columns() );

      const size_t jpos( prevMultiple( N, 2UL ) );
      BLAZE_INTERNAL_ASSERT( jpos <= N, "Invalid end calculation" );

      for( size_t i=0UL; i<M; ++i ) {
         for( size_t j=0UL; j<jpos; j+=2UL ) {
            A(i,j    ) = op( x[i], y[j  ] );
            A(i,j+1UL) = op( x[i], y[j+1] );
         }
         if( jpos < N ) {
            A(i,jpos) = op( x[i], y[jpos] );
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized assignment to row-major dense matrices*******************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized assignment of a dense vector-dense vector outer map expression to a
   //        row-major dense matrix.
   // \ingroup dense_matrix
   //
   // \param A The target left-hand side dense matrix.
   // \param x The left-hand side dense vector operand.
   // \param y The right-hand side dense vector operand.
   // \param op The custom binary operation.
   // \return void
   //
   // This function implements the vectorized assignment kernel for the dense vector-dense vector
   // outer map expression.
   */
   template< typename MT     // Type of the left-hand side target matrix
           , typename VT3    // Type of the left-hand side vector operand
           , typename VT4 >  // Type of the right-hand side vector operand
   static inline auto selectAssignKernel( MT& A, const VT3& x, const VT4& y, OP op )
      -> EnableIf_t< IsRowMajorMatrix_v<MT> && UseVectorizedKernel_v<MT,VT3,VT4> >
   {
      constexpr bool remainder( !IsSame_v<OP,Mult> || !IsPadded_v<MT> || !IsPadded_v<VT4> );

      const size_t M( A.rows() );
      const size_t N( A.columns() );

      const size_t jpos( remainder ? prevMultiple( N, SIMDSIZE ) : N );
      BLAZE_INTERNAL_ASSERT( jpos <= N, "Invalid end calculation" );

      auto xbegin( x.begin() );

      for( size_t i=0UL; i<M; ++i )
      {
         const auto x1( set( *xbegin ) );

         size_t j( 0UL );
         auto abegin( A.begin(i) );
         auto ybegin( y.begin()  );

         for( ; j<jpos; j+=SIMDSIZE, abegin+=SIMDSIZE, ybegin+=SIMDSIZE ) {
            abegin.store( op.load( x1, ybegin.load() ) );
         }
         for( ; remainder && j<N; ++j, ++abegin, ++ybegin ) {
            *abegin = op( *xbegin, *ybegin );
         }

         ++xbegin;
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to column-major dense matrices***************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense vector-dense vector outer map expression to a column-major
   //        dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side outer map expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense vector-dense
   // vector outer map expression to a column-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline void assign( DenseMatrix<MT,true>& lhs, const DVecDVecOuterExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT x( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense vector operand
      RT y( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (*lhs).rows()   , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == (*lhs).columns(), "Invalid vector size" );

      DVecDVecOuterExpr::selectAssignKernel( *lhs, x, y, rhs.op_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default assignment to column-major dense matrices*******************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default assignment of a dense vector-dense vector outer map expression to a
   //        column-major dense matrix.
   // \ingroup dense_matrix
   //
   // \param A The target left-hand side dense matrix.
   // \param x The left-hand side dense vector operand.
   // \param y The right-hand side dense vector operand.
   // \param op The custom binary operation.
   // \return void
   //
   // This function implements the default assignment kernel for the dense vector-dense vector
   // outer map expression.
   */
   template< typename MT     // Type of the left-hand side target matrix
           , typename VT3    // Type of the left-hand side vector operand
           , typename VT4 >  // Type of the right-hand side vector operand
   static inline auto selectAssignKernel( MT& A, const VT3& x, const VT4& y, OP op )
      -> EnableIf_t< IsColumnMajorMatrix_v<MT> && UseDefaultKernel_v<MT,VT3,VT4> >
   {
      const size_t M( A.rows() );
      const size_t N( A.columns() );

      const size_t ipos( prevMultiple( M, 2UL ) );
      BLAZE_INTERNAL_ASSERT( ipos <= M, "Invalid end calculation" );

      for( size_t j=0UL; j<N; ++j ) {
         for( size_t i=0UL; i<ipos; i+=2UL ) {
            A(i    ,j) = op( x[i  ], y[j] );
            A(i+1UL,j) = op( x[i+1], y[j] );
         }
         if( ipos < M ) {
            A(ipos,j) = op( x[ipos], y[j] );
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized assignment to column-major dense matrices****************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized assignment of a dense vector-dense vector outer map expression to a
   //        column-major dense matrix.
   // \ingroup dense_matrix
   //
   // \param A The target left-hand side dense matrix.
   // \param x The left-hand side dense vector operand.
   // \param y The right-hand side dense vector operand.
   // \param op The custom binary operation.
   // \return void
   //
   // This function implements the vectorized assignment kernel for the dense vector-dense vector
   // outer map expression.
   */
   template< typename MT     // Type of the left-hand side target matrix
           , typename VT3    // Type of the left-hand side vector operand
           , typename VT4 >  // Type of the right-hand side vector operand
   static inline auto selectAssignKernel( MT& A, const VT3& x, const VT4& y, OP op )
      -> EnableIf_t< IsColumnMajorMatrix_v<MT> && UseVectorizedKernel_v<MT,VT3,VT4> >
   {
      constexpr bool remainder( !IsSame_v<OP,Mult> || !IsPadded_v<MT> || !IsPadded_v<VT3> );

      const size_t M( A.rows() );
      const size_t N( A.columns() );

      const size_t ipos( remainder ? prevMultiple( M, SIMDSIZE ) : M );
      BLAZE_INTERNAL_ASSERT( ipos <= M, "Invalid end calculation" );

      auto ybegin( y.begin() );

      for( size_t j=0UL; j<N; ++j )
      {
         const auto y1( set( *ybegin ) );

         size_t i( 0UL );
         auto abegin( A.begin(j) );
         auto xbegin( x.begin()  );

         for( ; i<ipos; i+=SIMDSIZE, abegin+=SIMDSIZE, xbegin+=SIMDSIZE ) {
            abegin.store( op.load( xbegin.load(), y1 ) );
         }
         for( ; remainder && i<M; ++i, ++abegin, ++xbegin ) {
            *abegin = op( *xbegin, *ybegin );
         }

         ++ybegin;
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse matrices***************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense vector-dense vector outer map expression to a sparse matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side outer map expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense vector-dense
   // vector outer map expression to a sparse matrix.
   */
   template< typename MT  // Type of the target sparse matrix
           , bool SO >    // Storage order of the target sparse matrix
   friend inline void assign( SparseMatrix<MT,SO>& lhs, const DVecDVecOuterExpr& rhs )
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

   //**Addition assignment to row-major dense matrices*********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a dense vector-dense vector outer map to a row-major dense
   //        matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side outer map expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense vector-
   // dense vector outer map expression to a row-major dense matrix. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler
   // in case either of the two operands requires an intermediate evaluation.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline auto addAssign( DenseMatrix<MT,false>& lhs, const DVecDVecOuterExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT x( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense vector operand
      RT y( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (*lhs).rows()   , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == (*lhs).columns(), "Invalid vector size" );

      DVecDVecOuterExpr::selectAddAssignKernel( *lhs, x, y, rhs.op_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default addition assignment to row-major dense matrices*************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default addition assignment of a dense vector-dense vector outer map expression to
   //        a row-major dense matrix.
   // \ingroup dense_matrix
   //
   // \param A The target left-hand side dense matrix.
   // \param x The left-hand side dense vector operand.
   // \param y The right-hand side dense vector operand.
   // \param op The custom binary operation.
   // \return void
   //
   // This function implements the default addition assignment kernel for the dense vector-dense
   // vector outer map expression.
   */
   template< typename MT     // Type of the left-hand side target matrix
           , typename VT3    // Type of the left-hand side vector operand
           , typename VT4 >  // Type of the right-hand side vector operand
   static inline auto selectAddAssignKernel( MT& A, const VT3& x, const VT4& y, OP op )
      -> EnableIf_t< IsRowMajorMatrix_v<MT> && UseDefaultKernel_v<MT,VT3,VT4> >
   {
      const size_t M( A.rows() );
      const size_t N( A.columns() );

      const size_t jpos( prevMultiple( N, 2UL ) );
      BLAZE_INTERNAL_ASSERT( jpos <= N, "Invalid end calculation" );

      for( size_t i=0UL; i<M; ++i ) {
         for( size_t j=0UL; j<jpos; j+=2UL ) {
            A(i,j    ) += op( x[i], y[j    ] );
            A(i,j+1UL) += op( x[i], y[j+1UL] );
         }
         if( jpos < N ) {
            A(i,jpos) += op( x[i], y[jpos] );
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized addition assignment to row-major dense matrices**********************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized addition assignment of a dense vector-dense vector outer map expression
   //        to a row-major dense matrix.
   // \ingroup dense_matrix
   //
   // \param A The target left-hand side dense matrix.
   // \param x The left-hand side dense vector operand.
   // \param y The right-hand side dense vector operand.
   // \param op The custom binary operation.
   // \return void
   //
   // This function implements the vectorized addition assignment kernel for the dense vector-
   // dense vector outer map expression.
   */
   template< typename MT     // Type of the left-hand side target matrix
           , typename VT3    // Type of the left-hand side vector operand
           , typename VT4 >  // Type of the right-hand side vector operand
   static inline auto selectAddAssignKernel( MT& A, const VT3& x, const VT4& y, OP op )
      -> EnableIf_t< IsRowMajorMatrix_v<MT> && UseVectorizedKernel_v<MT,VT3,VT4> >
   {
      constexpr bool remainder( !IsSame_v<OP,Mult> || !IsPadded_v<MT> || !IsPadded_v<VT4> );

      const size_t M( A.rows() );
      const size_t N( A.columns() );

      const size_t jpos( remainder ? prevMultiple( N, SIMDSIZE ) : N );
      BLAZE_INTERNAL_ASSERT( jpos <= N, "Invalid end calculation" );

      auto xbegin( x.begin() );

      for( size_t i=0UL; i<M; ++i )
      {
         const auto x1( set( *xbegin ) );

         size_t j( 0UL );
         auto abegin( A.begin(i) );
         auto ybegin( y.begin()  );

         for( ; j<jpos; j+=SIMDSIZE, abegin+=SIMDSIZE, ybegin+=SIMDSIZE ) {
            abegin.store( abegin.load() + op.load( x1, ybegin.load() ) );
         }
         for( ; remainder && j<N; ++j, ++abegin, ++ybegin ) {
            *abegin += op( *xbegin, *ybegin );
         }

         ++xbegin;
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to column-major dense matrices******************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a dense vector-dense vector outer map to a column-major
   //        dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side outer map expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense vector-
   // dense vector outer map expression to a column-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline void addAssign( DenseMatrix<MT,true>& lhs, const DVecDVecOuterExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT x( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense vector operand
      RT y( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (*lhs).rows()   , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == (*lhs).columns(), "Invalid vector size" );

      DVecDVecOuterExpr::selectAddAssignKernel( *lhs, x, y, rhs.op_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default addition assignment to column dense matrices****************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default addition assignment of a dense vector-dense vector outer map expression to
   //        a column-major dense matrix.
   // \ingroup dense_matrix
   //
   // \param A The target left-hand side dense matrix.
   // \param x The left-hand side dense vector operand.
   // \param y The right-hand side dense vector operand.
   // \param op The custom binary operation.
   // \return void
   //
   // This function implements the default addition assignment kernel for the dense vector-dense
   // vector outer map expression.
   */
   template< typename MT     // Type of the left-hand side target matrix
           , typename VT3    // Type of the left-hand side vector operand
           , typename VT4 >  // Type of the right-hand side vector operand
   static inline auto selectAddAssignKernel( MT& A, const VT3& x, const VT4& y, OP op )
      -> EnableIf_t< IsColumnMajorMatrix_v<MT> && UseDefaultKernel_v<MT,VT3,VT4> >
   {
      const size_t M( A.rows() );
      const size_t N( A.columns() );

      const size_t ipos( prevMultiple( M, 2UL ) );
      BLAZE_INTERNAL_ASSERT( ipos <= M, "Invalid end calculation" );

      for( size_t j=0UL; j<N; ++j ) {
         for( size_t i=0UL; i<ipos; i+=2UL ) {
            A(i    ,j) += op( x[i    ], y[j] );
            A(i+1UL,j) += op( x[i+1UL], y[j] );
         }
         if( ipos < M ) {
            A(ipos,j) += op( x[ipos], y[j] );
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized addition assignment to column-major dense matrices*******************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized addition assignment of a dense vector-dense vector outer map expression
   //        to a column-major dense matrix.
   // \ingroup dense_matrix
   //
   // \param A The target left-hand side dense matrix.
   // \param x The left-hand side dense vector operand.
   // \param y The right-hand side dense vector operand.
   // \param op The custom binary operation.
   // \return void
   //
   // This function implements the vectorized addition assignment kernel for the dense vector-
   // dense vector outer map expression.
   */
   template< typename MT     // Type of the left-hand side target matrix
           , typename VT3    // Type of the left-hand side vector operand
           , typename VT4 >  // Type of the right-hand side vector operand
   static inline auto selectAddAssignKernel( MT& A, const VT3& x, const VT4& y, OP op )
      -> EnableIf_t< IsColumnMajorMatrix_v<MT> && UseVectorizedKernel_v<MT,VT3,VT4> >
   {
      constexpr bool remainder( !IsSame_v<OP,Mult> || !IsPadded_v<MT> || !IsPadded_v<VT3> );

      const size_t M( A.rows() );
      const size_t N( A.columns() );

      const size_t ipos( remainder ? prevMultiple( M, SIMDSIZE ) : M );
      BLAZE_INTERNAL_ASSERT( ipos <= M, "Invalid end calculation" );

      auto ybegin( y.begin() );

      for( size_t j=0UL; j<N; ++j )
      {
         const auto y1( set( *ybegin ) );

         size_t i( 0UL );
         auto abegin( A.begin(j) );
         auto xbegin( x.begin()  );

         for( ; i<ipos; i+=SIMDSIZE, abegin+=SIMDSIZE, xbegin+=SIMDSIZE ) {
            abegin.store( abegin.load() + op.load( xbegin.load(), y1 ) );
         }
         for( ; remainder && i<M; ++i, ++abegin, ++xbegin ) {
            *abegin += op( *xbegin, *ybegin );
         }

         ++ybegin;
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse matrices******************************************************
   // No special implementation for the addition assignment to sparse matrices.
   //**********************************************************************************************

   //**Subtraction assignment to row-major dense matrices******************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a dense vector-dense vector outer map to a row-major dense
   //        matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side outer map expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a dense
   // vector-dense vector outer map expression to a row-major dense matrix. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler
   // in case either of the two operands requires an intermediate evaluation.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline auto subAssign( DenseMatrix<MT,false>& lhs, const DVecDVecOuterExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT x( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense vector operand
      RT y( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (*lhs).rows()   , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == (*lhs).columns(), "Invalid vector size" );

      DVecDVecOuterExpr::selectSubAssignKernel( *lhs, x, y, rhs.op_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default subtraction assignment to row-major dense matrices**********************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default subtraction assignment of a dense vector-dense vector outer map expression
   //        to a row-major dense matrix.
   // \ingroup dense_matrix
   //
   // \param A The target left-hand side dense matrix.
   // \param x The left-hand side dense vector operand.
   // \param y The right-hand side dense vector operand.
   // \param op The custom binary operation.
   // \return void
   //
   // This function implements the default subtraction assignment kernel for the dense vector-
   // dense vector outer map expression.
   */
   template< typename MT     // Type of the left-hand side target matrix
           , typename VT3    // Type of the left-hand side vector operand
           , typename VT4 >  // Type of the right-hand side vector operand
   static inline auto selectSubAssignKernel( MT& A, const VT3& x, const VT4& y, OP op )
      -> EnableIf_t< IsRowMajorMatrix_v<MT> && UseDefaultKernel_v<MT,VT3,VT4> >
   {
      const size_t M( A.rows() );
      const size_t N( A.columns() );

      const size_t jpos( prevMultiple( N, 2UL ) );
      BLAZE_INTERNAL_ASSERT( jpos <= N, "Invalid end calculation" );

      for( size_t i=0UL; i<M; ++i ) {
         for( size_t j=0UL; j<jpos; j+=2UL ) {
            A(i,j    ) -= op( x[i], y[j    ] );
            A(i,j+1UL) -= op( x[i], y[j+1UL] );
         }
         if( jpos < N ) {
            A(i,jpos) -= op( x[i], y[jpos] );
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized subtraction assignment to row-major dense matrices*******************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized subtraction assignment of a dense vector-dense vector outer map expression
   //        to a row-major dense matrix.
   // \ingroup dense_matrix
   //
   // \param A The target left-hand side dense matrix.
   // \param x The left-hand side dense vector operand.
   // \param y The right-hand side dense vector operand.
   // \param op The custom binary operation.
   // \return void
   //
   // This function implements the vectorized subtraction assignment kernel for the dense vector-
   // dense vector outer map expression.
   */
   template< typename MT     // Type of the left-hand side target matrix
           , typename VT3    // Type of the left-hand side vector operand
           , typename VT4 >  // Type of the right-hand side vector operand
   static inline auto selectSubAssignKernel( MT& A, const VT3& x, const VT4& y, OP op )
      -> EnableIf_t< IsRowMajorMatrix_v<MT> && UseVectorizedKernel_v<MT,VT3,VT4> >
   {
      constexpr bool remainder( !IsSame_v<OP,Mult> || !IsPadded_v<MT> || !IsPadded_v<VT4> );

      const size_t M( A.rows() );
      const size_t N( A.columns() );

      const size_t jpos( remainder ? prevMultiple( N, SIMDSIZE ) : N );
      BLAZE_INTERNAL_ASSERT( jpos <= N, "Invalid end calculation" );

      auto xbegin( x.begin() );

      for( size_t i=0UL; i<M; ++i )
      {
         const auto x1( set( *xbegin ) );

         size_t j( 0UL );
         auto abegin( A.begin(i) );
         auto ybegin( y.begin()  );

         for( ; j<jpos; j+=SIMDSIZE, abegin+=SIMDSIZE, ybegin+=SIMDSIZE ) {
            abegin.store( abegin.load() - op.load( x1, ybegin.load() ) );
         }
         for( ; remainder && j<N; ++j, ++abegin, ++ybegin ) {
            *abegin -= op( *xbegin, *ybegin );
         }

         ++xbegin;
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to column-major dense matrices***************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a dense vector-dense vector outer map to a column-major
   //        dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side outer map expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a dense vector-
   // dense vector outer map expression to a column-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline void subAssign( DenseMatrix<MT,true>& lhs, const DVecDVecOuterExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT x( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense vector operand
      RT y( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (*lhs).rows()   , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == (*lhs).columns(), "Invalid vector size" );

      DVecDVecOuterExpr::selectSubAssignKernel( *lhs, x, y, rhs.op_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default subtraction assignment to column dense matrices*************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default subtraction assignment of a dense vector-dense vector outer map expression
   //        to a column-major dense matrix.
   // \ingroup dense_matrix
   //
   // \param A The target left-hand side dense matrix.
   // \param x The left-hand side dense vector operand.
   // \param y The right-hand side dense vector operand.
   // \param op The custom binary operation.
   // \return void
   //
   // This function implements the default subtraction assignment kernel for the dense vector-
   // dense vector outer map expression.
   */
   template< typename MT     // Type of the left-hand side target matrix
           , typename VT3    // Type of the left-hand side vector operand
           , typename VT4 >  // Type of the right-hand side vector operand
   static inline auto selectSubAssignKernel( MT& A, const VT3& x, const VT4& y, OP op )
      -> EnableIf_t< IsColumnMajorMatrix_v<MT> && UseDefaultKernel_v<MT,VT3,VT4> >
   {
      const size_t M( A.rows() );
      const size_t N( A.columns() );

      const size_t ipos( prevMultiple( M, 2UL ) );
      BLAZE_INTERNAL_ASSERT( ipos <= M, "Invalid end calculation" );

      for( size_t j=0UL; j<N; ++j ) {
         for( size_t i=0UL; i<ipos; i+=2UL ) {
            A(i    ,j) -= op( x[i    ], y[j] );
            A(i+1UL,j) -= op( x[i+1UL], y[j] );
         }
         if( ipos < M ) {
            A(ipos,j) -= op( x[ipos], y[j] );
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized subtraction assignment to column-major dense matrices****************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized subtraction assignment of a dense vector-dense vector outer map expression
   //        to a column-major dense matrix.
   // \ingroup dense_matrix
   //
   // \param A The target left-hand side dense matrix.
   // \param x The left-hand side dense vector operand.
   // \param y The right-hand side dense vector operand.
   // \param op The custom binary operation.
   // \return void
   //
   // This function implements the vectorized subtraction assignment kernel for the dense vector-
   // dense vector outer map expression.
   */
   template< typename MT     // Type of the left-hand side target matrix
           , typename VT3    // Type of the left-hand side vector operand
           , typename VT4 >  // Type of the right-hand side vector operand
   static inline auto selectSubAssignKernel( MT& A, const VT3& x, const VT4& y, OP op )
      -> EnableIf_t< IsColumnMajorMatrix_v<MT> && UseVectorizedKernel_v<MT,VT3,VT4> >
   {
      constexpr bool remainder( !IsSame_v<OP,Mult> || !IsPadded_v<MT> || !IsPadded_v<VT3> );

      const size_t M( A.rows() );
      const size_t N( A.columns() );

      const size_t ipos( remainder ? prevMultiple( M, SIMDSIZE ) : M );
      BLAZE_INTERNAL_ASSERT( ipos <= M, "Invalid end calculation" );

      auto ybegin( y.begin() );

      for( size_t j=0UL; j<N; ++j )
      {
         const auto y1( set( *ybegin ) );

         size_t i( 0UL );
         auto abegin( A.begin(j) );
         auto xbegin( x.begin()  );

         for( ; i<ipos; i+=SIMDSIZE, abegin+=SIMDSIZE, xbegin+=SIMDSIZE ) {
            abegin.store( abegin.load() - op.load( xbegin.load(), y1 ) );
         }
         for( ; remainder && i<M; ++i, ++abegin, ++xbegin ) {
            *abegin -= op( *xbegin, *ybegin );
         }

         ++ybegin;
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse matrices***************************************************
   // No special implementation for the subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**Schur product assignment to row-major dense matrices****************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Schur product assignment of a dense vector-dense vector outer map to a row-major
   //        dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side outer map expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a dense
   // vector-dense vector outer map expression to a row-major dense matrix. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler
   // in case either of the two operands requires an intermediate evaluation.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline auto schurAssign( DenseMatrix<MT,false>& lhs, const DVecDVecOuterExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT x( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense vector operand
      RT y( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (*lhs).rows()   , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == (*lhs).columns(), "Invalid vector size" );

      DVecDVecOuterExpr::selectSchurAssignKernel( *lhs, x, y, rhs.op_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default Schur product assignment to row-major dense matrices********************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default Schur product assignment of a dense vector-dense vector outer map expression
   //        to a row-major dense matrix.
   // \ingroup dense_matrix
   //
   // \param A The target left-hand side dense matrix.
   // \param x The left-hand side dense vector operand.
   // \param y The right-hand side dense vector operand.
   // \param op The custom binary operation.
   // \return void
   //
   // This function implements the default Schur product assignment kernel for the dense vector-
   // dense vector outer map expression.
   */
   template< typename MT     // Type of the left-hand side target matrix
           , typename VT3    // Type of the left-hand side vector operand
           , typename VT4 >  // Type of the right-hand side vector operand
   static inline auto selectSchurAssignKernel( MT& A, const VT3& x, const VT4& y, OP op )
      -> EnableIf_t< IsRowMajorMatrix_v<MT> && UseDefaultKernel_v<MT,VT3,VT4> >
   {
      const size_t M( A.rows() );
      const size_t N( A.columns() );

      const size_t jpos( prevMultiple( N, 2UL ) );
      BLAZE_INTERNAL_ASSERT( jpos <= N, "Invalid end calculation" );

      for( size_t i=0UL; i<M; ++i ) {
         for( size_t j=0UL; j<jpos; j+=2UL ) {
            A(i,j    ) *= op( x[i], y[j    ] );
            A(i,j+1UL) *= op( x[i], y[j+1UL] );
         }
         if( jpos < N ) {
            A(i,jpos) *= op( x[i], y[jpos] );
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized Schur product assignment to row-major dense matrices*****************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized Schur product assignment of a dense vector-dense vector outer map
   //        expression to a row-major dense matrix.
   // \ingroup dense_matrix
   //
   // \param A The target left-hand side dense matrix.
   // \param x The left-hand side dense vector operand.
   // \param y The right-hand side dense vector operand.
   // \param op The custom binary operation.
   // \return void
   //
   // This function implements the vectorized Schur product assignment kernel for the dense vector-
   // dense vector outer map expression.
   */
   template< typename MT     // Type of the left-hand side target matrix
           , typename VT3    // Type of the left-hand side vector operand
           , typename VT4 >  // Type of the right-hand side vector operand
   static inline auto selectSchurAssignKernel( MT& A, const VT3& x, const VT4& y, OP op )
      -> EnableIf_t< IsRowMajorMatrix_v<MT> && UseVectorizedKernel_v<MT,VT3,VT4> >
   {
      constexpr bool remainder( !IsSame_v<OP,Mult> || !IsPadded_v<MT> || !IsPadded_v<VT4> );

      const size_t M( A.rows() );
      const size_t N( A.columns() );

      const size_t jpos( remainder ? prevMultiple( N, SIMDSIZE ) : N );
      BLAZE_INTERNAL_ASSERT( jpos <= N, "Invalid end calculation" );

      auto xbegin( x.begin() );

      for( size_t i=0UL; i<M; ++i )
      {
         const auto x1( set( *xbegin ) );

         size_t j( 0UL );
         auto abegin( A.begin(i) );
         auto ybegin( y.begin()  );

         for( ; j<jpos; j+=SIMDSIZE, abegin+=SIMDSIZE, ybegin+=SIMDSIZE ) {
            abegin.store( abegin.load() * op.load( x1, ybegin.load() ) );
         }
         for( ; remainder && j<N; ++j, ++abegin, ++ybegin ) {
            *abegin *= op( *xbegin, *ybegin );
         }

         ++xbegin;
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Schur product assignment to column-major dense matrices*************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Schur product assignment of a dense vector-dense vector outer map to a column-major
   //        dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side outer map expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a dense
   // vector-dense vector outer map expression to a column-major dense matrix.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline void schurAssign( DenseMatrix<MT,true>& lhs, const DVecDVecOuterExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT x( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense vector operand
      RT y( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (*lhs).rows()   , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == (*lhs).columns(), "Invalid vector size" );

      DVecDVecOuterExpr::selectSchurAssignKernel( *lhs, x, y, rhs.op_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default Schur product assignment to column dense matrices***********************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default Schur product assignment of a dense vector-dense vector outer map expression
   //        to a column-major dense matrix.
   // \ingroup dense_matrix
   //
   // \param A The target left-hand side dense matrix.
   // \param x The left-hand side dense vector operand.
   // \param y The right-hand side dense vector operand.
   // \param op The custom binary operation.
   // \return void
   //
   // This function implements the default Schur product assignment kernel for the dense vector-
   // dense vector outer map expression.
   */
   template< typename MT     // Type of the left-hand side target matrix
           , typename VT3    // Type of the left-hand side vector operand
           , typename VT4 >  // Type of the right-hand side vector operand
   static inline auto selectSchurAssignKernel( MT& A, const VT3& x, const VT4& y, OP op )
      -> EnableIf_t< IsColumnMajorMatrix_v<MT> && UseDefaultKernel_v<MT,VT3,VT4> >
   {
      const size_t M( A.rows() );
      const size_t N( A.columns() );

      const size_t ipos( prevMultiple( M, 2UL ) );
      BLAZE_INTERNAL_ASSERT( ipos <= M, "Invalid end calculation" );

      for( size_t j=0UL; j<N; ++j ) {
         for( size_t i=0UL; i<ipos; i+=2UL ) {
            A(i    ,j) *= op( x[i    ], y[j] );
            A(i+1UL,j) *= op( x[i+1UL], y[j] );
         }
         if( ipos < M ) {
            A(ipos,j) *= op( x[ipos], y[j] );
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized Schur product assignment to column-major dense matrices**************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized Schur product assignment of a dense vector-dense vector outer map
   //        expression to a column-major dense matrix.
   // \ingroup dense_matrix
   //
   // \param A The target left-hand side dense matrix.
   // \param x The left-hand side dense vector operand.
   // \param y The right-hand side dense vector operand.
   // \param op The custom binary operation.
   // \return void
   //
   // This function implements the vectorized Schur product assignment kernel for the dense vector-
   // dense vector outer map expression.
   */
   template< typename MT     // Type of the left-hand side target matrix
           , typename VT3    // Type of the left-hand side vector operand
           , typename VT4 >  // Type of the right-hand side vector operand
   static inline auto selectSchurAssignKernel( MT& A, const VT3& x, const VT4& y, OP op )
      -> EnableIf_t< IsColumnMajorMatrix_v<MT> && UseVectorizedKernel_v<MT,VT3,VT4> >
   {
      constexpr bool remainder( !IsSame_v<OP,Mult> || !IsPadded_v<MT> || !IsPadded_v<VT3> );

      const size_t M( A.rows() );
      const size_t N( A.columns() );

      const size_t ipos( remainder ? prevMultiple( M, SIMDSIZE ) : M );
      BLAZE_INTERNAL_ASSERT( ipos <= M, "Invalid end calculation" );

      auto ybegin( y.begin() );

      for( size_t j=0UL; j<N; ++j )
      {
         const auto y1( set( *ybegin ) );

         size_t i( 0UL );
         auto abegin( A.begin(j) );
         auto xbegin( x.begin()  );

         for( ; i<ipos; i+=SIMDSIZE, abegin+=SIMDSIZE, xbegin+=SIMDSIZE ) {
            abegin.store( abegin.load() * op.load( xbegin.load(), y1 ) );
         }
         for( ; remainder && i<M; ++i, ++abegin, ++xbegin ) {
            *abegin *= op( *xbegin, *ybegin );
         }

         ++ybegin;
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
   /*!\brief SMP assignment of a dense vector-dense vector outer map expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side outer map expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a dense vector-dense
   // vector outer map expression to a dense matrix. Due to the explicit application of the
   // SFINAE principle, this function can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto smpAssign( DenseMatrix<MT,SO>& lhs, const DVecDVecOuterExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT x( rhs.lhs_ );  // Evaluation of the left-hand side dense vector operand
      RT y( rhs.rhs_ );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (*lhs).rows()   , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == (*lhs).columns(), "Invalid vector size" );

      smpAssign( *lhs, map( x, y, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to sparse matrices***********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a dense vector-dense vector outer map exression to a sparse matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side outer map expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a dense vector-dense
   // vector outer map expression to a sparse matrix. Due to the explicit application of the
   // SFINAE principle, this function can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename MT  // Type of the target sparse matrix
           , bool SO >    // Storage order of the target sparse matrix
   friend inline auto smpAssign( SparseMatrix<MT,SO>& lhs, const DVecDVecOuterExpr& rhs )
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
   /*!\brief SMP addition assignment of a dense vector-dense vector outer map expression to a
   //        dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side outer map expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a dense
   // vector-dense vector outer map expression to a dense matrix. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler
   // in case the expression specific parallel evaluation strategy is selected.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline auto smpAddAssign( DenseMatrix<MT,false>& lhs, const DVecDVecOuterExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT x( rhs.lhs_ );  // Evaluation of the left-hand side dense vector operand
      RT y( rhs.rhs_ );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (*lhs).rows()   , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == (*lhs).columns(), "Invalid vector size" );

      smpAddAssign( *lhs, map( x, y, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to sparse matrices**************************************************
   // No special implementation for the SMP addition assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP subtraction assignment to dense matrices************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a dense vector-dense vector outer map expression to a
   //        dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side outer map expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a dense
   // vector-dense vector outer map expression to a dense matrix. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler
   // in case the expression specific parallel evaluation strategy is selected.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline auto smpSubAssign( DenseMatrix<MT,false>& lhs, const DVecDVecOuterExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT x( rhs.lhs_ );  // Evaluation of the left-hand side dense vector operand
      RT y( rhs.rhs_ );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (*lhs).rows()   , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == (*lhs).columns(), "Invalid vector size" );

      smpSubAssign( *lhs, map( x, y, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to sparse matrices***********************************************
   // No special implementation for the SMP subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP Schur product assignment to dense matrices**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP Schur product assignment of a dense vector-dense vector outer map expression to
   //        a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side outer map expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized SMP Schur product assignment of a dense
   // vector-dense vector outer map expression to a dense matrix. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler
   // in case the expression specific parallel evaluation strategy is selected.
   */
   template< typename MT >  // Type of the target dense matrix
   friend inline auto smpSchurAssign( DenseMatrix<MT,false>& lhs, const DVecDVecOuterExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT x( rhs.lhs_ );  // Evaluation of the left-hand side dense vector operand
      RT y( rhs.rhs_ );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size() , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (*lhs).rows()   , "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == (*lhs).columns(), "Invalid vector size" );

      smpSchurAssign( *lhs, map( x, y, rhs.op_ ) );
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
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( VT1 );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( VT2 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( VT1 );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE   ( VT2 );
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
/*!\brief Pairwise (outer) evaluation of the given binary operation on the elements of the dense
//        vectors \a lhs and \a rhs.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense vector for the outer map.
// \param rhs The right-hand side transpose dense vector for the outer map.
// \param op The custom, binary operation.
// \return The binary operation applied to each pair of elements of \a lhs and \a rhs.
//
// The \a map() function evaluates the given binary operation on each pair of elements of the
// input vectors \a lhs and \a rhs. The function returns an expression representing this
// operation.\n
// The following example demonstrates the use of the \a map() function:

   \code
   using blaze::columnVector;
   using blaze::rowMajor;

   blaze::DynamicVector<double,columnVector> a, b;
   blaze::DynamicMatrix<double,rowMajor> A;
   // ... Resizing and initialization
   A = map( a, trans(b), []( double x, double y ){ return a + b; );
   \endcode
*/
template< typename VT1   // Type of the left-hand side dense vector
        , typename VT2   // Type of the right-hand side dense vector
        , typename OP >  // Type of the custom operation
inline decltype(auto)
   map( const DenseVector<VT1,false>& lhs, const DenseVector<VT2,true>& rhs, OP op )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const DVecDVecOuterExpr<VT1,VT2,OP>;
   return ReturnType( *lhs, *rhs, std::move(op) );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Addition operator for the outer sum of two dense vectors (\f$ A=\vec{b}+\vec{c}^T \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense vector for the outer sum.
// \param rhs The right-hand side transpose dense vector for the outer sum.
// \return The resulting dense matrix.
//
// This operator represents the outer sum between a dense vector and a transpose dense
// vector:

   \code
   blaze::DynamicVector<double> a, b, c;
   // ... Resizing and initialization
   c = a + trans(b);
   \endcode

// The operator returns an expression representing a dense matrix of the higher-order element
// type of the two involved element types \a VT1::ElementType and \a VT2::ElementType.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
inline decltype(auto)
   operator+( const DenseVector<VT1,false>& lhs, const DenseVector<VT2,true>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return map( *lhs, *rhs, Add() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition operator for the outer difference of two dense vectors
//        (\f$ A=\vec{b}-\vec{c}^T \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense vector for the outer difference.
// \param rhs The right-hand side transpose dense vector for the outer difference.
// \return The resulting dense matrix.
//
// This operator represents the outer difference between a dense vector and a transpose dense
// vector:

   \code
   blaze::DynamicVector<double> a, b, c;
   // ... Resizing and initialization
   c = a - trans(b);
   \endcode

// The operator returns an expression representing a dense matrix of the higher-order element
// type of the two involved element types \a VT1::ElementType and \a VT2::ElementType.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
inline decltype(auto)
   operator-( const DenseVector<VT1,false>& lhs, const DenseVector<VT2,true>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return map( *lhs, *rhs, Sub() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the outer product of two dense vectors
//        (\f$ A=\vec{b}*\vec{c}^T \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense vector for the outer product.
// \param rhs The right-hand side transpose dense vector for the outer product.
// \return The resulting dense matrix.
//
// This operator represents the outer product between a dense vector and a transpose dense
// vector:

   \code
   using blaze::columnVector;
   using blaze::rowMajor;

   blaze::DynamicVector<double,columnVector> a, b;
   blaze::DynamicMatrix<double,rowMajor> A;
   // ... Resizing and initialization
   A = a * trans(b);
   \endcode

// The operator returns an expression representing a dense matrix of the higher-order element
// type of the two involved element types \a VT1::ElementType and \a VT2::ElementType.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
inline decltype(auto)
   operator*( const DenseVector<VT1,false>& lhs, const DenseVector<VT2,true>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return map( *lhs, *rhs, Mult() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division operator for the outer quotient of two dense vectors
//        (\f$ A=\vec{b}/\vec{c}^T \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense vector for the outer quotient.
// \param rhs The right-hand side transpose dense vector for the outer quotient.
// \return The resulting dense matrix.
//
// This operator represents the outer quotient between a dense vector and a transpose dense
// vector:

   \code
   using blaze::columnVector;
   using blaze::rowMajor;

   blaze::DynamicVector<double,columnVector> a, b;
   blaze::DynamicMatrix<double,rowMajor> A;
   // ... Resizing and initialization
   A = a / trans(b);
   \endcode

// The operator returns an expression representing a dense matrix of the higher-order element
// type of the two involved element types \a VT1::ElementType and \a VT2::ElementType.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
inline decltype(auto)
   operator/( const DenseVector<VT1,false>& lhs, const DenseVector<VT2,true>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return map( *lhs, *rhs, Div() );
}
//*************************************************************************************************




//=================================================================================================
//
//  SIZE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename VT2, typename OP >
struct Size< DVecDVecOuterExpr<VT1,VT2,OP>, 0UL >
   : public Size<VT1,0UL>
{};

template< typename VT1, typename VT2, typename OP >
struct Size< DVecDVecOuterExpr<VT1,VT2,OP>, 1UL >
   : public Size<VT2,0UL>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISALIGNED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename VT2, typename OP >
struct IsAligned< DVecDVecOuterExpr<VT1,VT2,OP> >
   : public BoolConstant< IsAligned_v<VT1> && IsAligned_v<VT2> >
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
template< typename VT1, typename VT2 >
struct IsPadded< DVecDVecOuterExpr<VT1,VT2,Mult> >
   : public IsPadded<VT2>
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
