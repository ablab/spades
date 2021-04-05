//=================================================================================================
/*!
//  \file blaze/math/expressions/DMatDMatMapExpr.h
//  \brief Header file for the dense matrix/dense matrix map expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DMATDMATMAPEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DMATDMATMAPEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <utility>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/MatMatMapExpr.h>
#include <blaze/math/functors/And.h>
#include <blaze/math/functors/Atan2.h>
#include <blaze/math/functors/Bitand.h>
#include <blaze/math/functors/Bitor.h>
#include <blaze/math/functors/Bitxor.h>
#include <blaze/math/functors/Hypot.h>
#include <blaze/math/functors/Join.h>
#include <blaze/math/functors/MakePair.h>
#include <blaze/math/functors/Max.h>
#include <blaze/math/functors/Min.h>
#include <blaze/math/functors/Or.h>
#include <blaze/math/functors/Pow.h>
#include <blaze/math/functors/ShiftLV.h>
#include <blaze/math/functors/ShiftRV.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/SIMD.h>
#include <blaze/math/traits/MapTrait.h>
#include <blaze/math/typetraits/HasLoad.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/IsPadded.h>
#include <blaze/math/typetraits/IsPaddingEnabled.h>
#include <blaze/math/typetraits/IsSIMDEnabled.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/system/HostDevice.h>
#include <blaze/system/Inline.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/HasMember.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DMATDMATMAPEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for the dense matrix-dense matrix map() function.
// \ingroup dense_matrix_expression
//
// The DMatDMatMapExpr class represents the compile time expression for the pairwise evaluation
// of a binary custom operation on the elements of two dense matrices with identical storage order
// via the map() function.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , typename OP   // Type of the custom operation
        , bool SO >     // Storage order
class DMatDMatMapExpr
   : public MatMatMapExpr< DenseMatrix< DMatDMatMapExpr<MT1,MT2,OP,SO>, SO > >
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   using RT1 = ResultType_t<MT1>;     //!< Result type of the left-hand side dense matrix expression.
   using RT2 = ResultType_t<MT2>;     //!< Result type of the right-hand side dense matrix expression.
   using ET1 = ElementType_t<MT1>;    //!< Element type of the left-hand side dense matrix expression.
   using ET2 = ElementType_t<MT2>;    //!< Element type of the right-hand side dense matrix expression.
   using RN1 = ReturnType_t<MT1>;     //!< Return type of the left-hand side dense matrix expression.
   using RN2 = ReturnType_t<MT2>;     //!< Return type of the right-hand side dense matrix expression.
   using CT1 = CompositeType_t<MT1>;  //!< Composite type of the left-hand side dense matrix expression.
   using CT2 = CompositeType_t<MT2>;  //!< Composite type of the right-hand side dense matrix expression.
   //**********************************************************************************************

   //**Serial evaluation strategy******************************************************************
   //! Compilation switch for the serial evaluation strategy of the map expression.
   /*! The \a useAssign compile time constant expression represents a compilation switch for
       the serial evaluation strategy of the map expression. In case either of the two dense
       matrix operands requires an intermediate evaluation, \a useAssign will be set to \a true
       and the addition expression will be evaluated via the \a assign function family. Otherwise
       \a useAssign will be set to \a false and the expression will be evaluated via the subscript
       operator. */
   static constexpr bool useAssign = ( RequiresEvaluation_v<MT1> || RequiresEvaluation_v<MT2> );

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
       In case at least one of the two dense matrix operands is not SMP assignable and at least
       one of the two operands requires an intermediate evaluation, the variable is set to \a true
       and the expression specific evaluation strategy is selected. Otherwise the variable is set
       to \a false and the default strategy is chosen. */
   template< typename MT >
   static constexpr bool UseSMPAssign_v =
      ( ( !MT1::smpAssignable || !MT2::smpAssignable ) && useAssign );
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this DMatDMatMapExpr instance.
   using This = DMatDMatMapExpr<MT1,MT2,OP,SO>;

   //! Base type of this DMatDMatMapExpr instance.
   using BaseType = MatMatMapExpr< DenseMatrix<This,SO> >;

   using ResultType    = MapTrait_t<RT1,RT2,OP>;       //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<ResultType>;   //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<ResultType>;    //!< Resulting element type.

   //! Return type for expression template evaluations.
   using ReturnType = decltype( std::declval<OP>()( std::declval<RN1>(), std::declval<RN2>() ) );

   //! Data type for composite expression templates.
   using CompositeType = If_t< useAssign, const ResultType, const DMatDMatMapExpr& >;

   //! Composite type of the left-hand side dense matrix expression.
   using LeftOperand = If_t< IsExpression_v<MT1>, const MT1, const MT1& >;

   //! Composite type of the right-hand side dense matrix expression.
   using RightOperand = If_t< IsExpression_v<MT2>, const MT2, const MT2& >;

   //! Data type of the custom binary operation.
   using Operation = OP;

   //! Type for the assignment of the left-hand side dense matrix operand.
   using LT = If_t< RequiresEvaluation_v<MT1>, const RT1, CT1 >;

   //! Type for the assignment of the right-hand side dense matrix operand.
   using RT = If_t< RequiresEvaluation_v<MT2>, const RT2, CT2 >;
   //**********************************************************************************************

   //**ConstIterator class definition**************************************************************
   /*!\brief Iterator over the elements of the dense matrix map expression.
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

      //! ConstIterator type of the left-hand side dense matrix expression.
      using LeftIteratorType = ConstIterator_t<MT1>;

      //! ConstIterator type of the right-hand side dense matrix expression.
      using RightIteratorType = ConstIterator_t<MT2>;
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
         left_  += inc;
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
         left_  -= dec;
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
         ++left_;
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
         return ConstIterator( left_++, right_++, op_ );
      }
      //*******************************************************************************************

      //**Prefix decrement operator****************************************************************
      /*!\brief Pre-decrement operator.
      //
      // \return Reference to the decremented iterator.
      */
      inline BLAZE_DEVICE_CALLABLE ConstIterator& operator--() {
         --left_;
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
         return ConstIterator( left_--, right_--, op_ );
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
      /*!\brief Access to the SIMD elements of the matrix.
      //
      // \return The resulting SIMD element.
      */
      inline auto load() const noexcept {
         return op_.load( left_.load(), right_.load() );
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator==( const ConstIterator& rhs ) const {
         return left_ == rhs.left_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator!=( const ConstIterator& rhs ) const {
         return left_ != rhs.left_;
      }
      //*******************************************************************************************

      //**Less-than operator***********************************************************************
      /*!\brief Less-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller, \a false if not.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator<( const ConstIterator& rhs ) const {
         return left_ < rhs.left_;
      }
      //*******************************************************************************************

      //**Greater-than operator********************************************************************
      /*!\brief Greater-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater, \a false if not.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator>( const ConstIterator& rhs ) const {
         return left_ > rhs.left_;
      }
      //*******************************************************************************************

      //**Less-or-equal-than operator**************************************************************
      /*!\brief Less-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller or equal, \a false if not.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator<=( const ConstIterator& rhs ) const {
         return left_ <= rhs.left_;
      }
      //*******************************************************************************************

      //**Greater-or-equal-than operator***********************************************************
      /*!\brief Greater-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater or equal, \a false if not.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator>=( const ConstIterator& rhs ) const {
         return left_ >= rhs.left_;
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two iterators.
      //
      // \param rhs The right-hand side iterator.
      // \return The number of elements between the two iterators.
      */
      inline BLAZE_DEVICE_CALLABLE DifferenceType operator-( const ConstIterator& rhs ) const {
         return left_ - rhs.left_;
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
         return ConstIterator( it.left_ + inc, it.right_ + inc, it.op_ );
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
         return ConstIterator( it.left_ - dec, it.right_ - dec, it.op_ );
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
      ( MT1::simdEnabled && MT2::simdEnabled &&
        If_t< HasSIMDEnabled_v<OP>, GetSIMDEnabled<OP,ET1,ET2>, HasLoad<OP> >::value );

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = ( MT1::smpAssignable && MT2::smpAssignable );
   //**********************************************************************************************

   //**SIMD properties*****************************************************************************
   //! The number of elements packed within a single SIMD element.
   static constexpr size_t SIMDSIZE = SIMDTrait<ElementType>::size;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DMatDMatMapExpr class.
   //
   // \param lhs The left-hand side dense matrix operand of the map expression.
   // \param rhs The right-hand side dense matrix operand of the map expression.
   // \param op The custom binary operation.
   */
   inline DMatDMatMapExpr( const MT1& lhs, const MT2& rhs, OP op ) noexcept
      : lhs_( lhs )            // Left-hand side dense matrix of the map expression
      , rhs_( rhs )            // Right-hand side dense matrix of the map expression
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
      BLAZE_INTERNAL_ASSERT( i < lhs_.rows()   , "Invalid row access index"    );
      BLAZE_INTERNAL_ASSERT( j < lhs_.columns(), "Invalid column access index" );
      return op_( lhs_(i,j), rhs_(i,j) );
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

   //**Load function*******************************************************************************
   /*!\brief Access to the SIMD elements of the matrix.
   //
   // \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
   // \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
   // \return Reference to the accessed values.
   */
   BLAZE_ALWAYS_INLINE auto load( size_t i, size_t j ) const noexcept {
      BLAZE_INTERNAL_ASSERT( i < lhs_.rows()   , "Invalid row access index"    );
      BLAZE_INTERNAL_ASSERT( j < lhs_.columns(), "Invalid column access index" );
      BLAZE_INTERNAL_ASSERT( !SO || ( i % SIMDSIZE == 0UL ), "Invalid row access index"    );
      BLAZE_INTERNAL_ASSERT( SO  || ( j % SIMDSIZE == 0UL ), "Invalid column access index" );
      return op_.load( lhs_.load(i,j), rhs_.load(i,j) );
   }
   //**********************************************************************************************

   //**Begin function******************************************************************************
   /*!\brief Returns an iterator to the first non-zero element of row/column \a i.
   //
   // \param i The row/column index.
   // \return Iterator to the first non-zero element of row/column \a i.
   */
   inline ConstIterator begin( size_t i ) const {
      return ConstIterator( lhs_.begin(i), rhs_.begin(i), op_ );
   }
   //**********************************************************************************************

   //**End function********************************************************************************
   /*!\brief Returns an iterator just past the last non-zero element of row/column \a i.
   //
   // \param i The row/column index.
   // \return Iterator just past the last non-zero element of row/column \a i.
   */
   inline ConstIterator end( size_t i ) const {
      return ConstIterator( lhs_.end(i), rhs_.end(i), op_ );
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
   /*!\brief Returns the right-hand side dense matrix operand.
   //
   // \return The right-hand side dense matrix operand.
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
      return ( IsExpression_v<MT1> && lhs_.canAlias( alias ) ) ||
             ( IsExpression_v<MT2> && rhs_.canAlias( alias ) );
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
      return lhs_.canSMPAssign() && rhs_.canSMPAssign();
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  lhs_;  //!< Left-hand side dense matrix of the map expression.
   RightOperand rhs_;  //!< Right-hand side dense matrix of the map expression.
   Operation    op_;   //!< The custom binary operation.
   //**********************************************************************************************

   //**Assignment to dense matrices****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense matrix-dense matrix map expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix-dense
   // matrix map expression to a dense matrix. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case either of the two
   // operands requires an intermediate evaluation.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO2 >   // Storage order of the target dense matrix
   friend inline auto assign( DenseMatrix<MT,SO2>& lhs, const DMatDMatMapExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense matrix operand
      RT B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == (*lhs).columns()  , "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (*lhs).columns()  , "Invalid number of columns" );

      assign( *lhs, map( A, B, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse matrices***************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense matrix-dense matrix map expression to a sparse matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side map expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix-dense
   // matrix map expression to a sparse matrix. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case either of the two
   // operands requires an intermediate evaluation.
   */
   template< typename MT  // Type of the target sparse matrix
           , bool SO2 >   // Storage order of the target sparse matrix
   friend inline auto assign( SparseMatrix<MT,SO2>& lhs, const DMatDMatMapExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT> >
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

   //**Addition assignment to dense matrices*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a dense matrix-dense matrix map expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense
   // matrix-dense matrix map expression to a dense matrix. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case
   // either of the two operands requires an intermediate evaluation.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO2 >   // Storage order of the target dense matrix
   friend inline auto addAssign( DenseMatrix<MT,SO2>& lhs, const DMatDMatMapExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense matrix operand
      RT B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == (*lhs).columns()  , "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (*lhs).columns()  , "Invalid number of columns" );

      addAssign( *lhs, map( A, B, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse matrices******************************************************
   // No special implementation for the addition assignment to sparse matrices.
   //**********************************************************************************************

   //**Subtraction assignment to dense matrices****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a dense matrix-dense matrix map expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a dense
   // matrix-dense matrix map expression to a dense matrix. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case
   // either of the two operands requires an intermediate evaluation.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO2 >   // Storage order of the target dense matrix
   friend inline auto subAssign( DenseMatrix<MT,SO2>& lhs, const DMatDMatMapExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense matrix operand
      RT B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == (*lhs).columns()  , "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (*lhs).columns()  , "Invalid number of columns" );

      subAssign( *lhs, map( A, B, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse matrices***************************************************
   // No special implementation for the subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**Schur product assignment to dense matrices**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Schur product assignment of a dense matrix-dense matrix map expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a dense
   // matrix-dense matrix map expression to a dense matrix. Due to the explicit application of
   // the SFINAE principle, this function can only be selected by the compiler in case either
   // of the two operands requires an intermediate evaluation.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO2 >   // Storage order of the target dense matrix
   friend inline auto schurAssign( DenseMatrix<MT,SO2>& lhs, const DMatDMatMapExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense matrix operand
      RT B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == (*lhs).columns()  , "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (*lhs).columns()  , "Invalid number of columns" );

      schurAssign( *lhs, map( A, B, rhs.op_ ) );
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
   /*!\brief SMP assignment of a dense matrix-dense matrix map expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a dense matrix-dense
   // matrix map expression to a dense matrix. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO2 >   // Storage order of the target dense matrix
   friend inline auto smpAssign( DenseMatrix<MT,SO2>& lhs, const DMatDMatMapExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( rhs.lhs_ );  // Evaluation of the left-hand side dense matrix operand
      RT B( rhs.rhs_ );  // Evaluation of the right-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == (*lhs).columns()  , "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (*lhs).columns()  , "Invalid number of columns" );

      smpAssign( *lhs, map( A, B, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to sparse matrices***********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a dense matrix-dense matrix map expression to a sparse matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side map expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a dense matrix-dense
   // matrix map expression to a sparse matrix. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename MT  // Type of the target sparse matrix
           , bool SO2 >   // Storage order of the target sparse matrix
   friend inline auto smpAssign( SparseMatrix<MT,SO2>& lhs, const DMatDMatMapExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT> >
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

      const TmpType tmp( rhs );
      smpAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to dense matrices***************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a dense matrix-dense matrix map expression to a dense
   //        matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a dense
   // matrix-dense matrix map expression to a dense matrix. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO2 >   // Storage order of the target dense matrix
   friend inline auto smpAddAssign( DenseMatrix<MT,SO2>& lhs, const DMatDMatMapExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( rhs.lhs_ );  // Evaluation of the left-hand side dense matrix operand
      RT B( rhs.rhs_ );  // Evaluation of the right-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == (*lhs).columns()  , "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (*lhs).columns()  , "Invalid number of columns" );

      smpAddAssign( *lhs, map( A, B, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to sparse matrices**************************************************
   // No special implementation for the SMP addition assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP subtraction assignment to dense matrices************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a dense matrix-dense matrix map expression to a dense
   //        matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a dense
   // matrix-dense matrix map expression to a dense matrix. Due to the explicit application of
   // the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO2 >   // Storage order of the target dense matrix
   friend inline auto smpSubAssign( DenseMatrix<MT,SO2>& lhs, const DMatDMatMapExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( rhs.lhs_ );  // Evaluation of the left-hand side dense matrix operand
      RT B( rhs.rhs_ );  // Evaluation of the right-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == (*lhs).columns()  , "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (*lhs).columns()  , "Invalid number of columns" );

      smpSubAssign( *lhs, map( A, B, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to sparse matrices***********************************************
   // No special implementation for the SMP subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP Schur product assignment to dense matrices**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP Schur product assignment of a dense matrix-dense matrix map expression to a
   //        dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized SMP Schur product assignment of a
   // dense matrix-dense matrix map expression to a dense matrix. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO2 >   // Storage order of the target dense matrix
   friend inline auto smpSchurAssign( DenseMatrix<MT,SO2>& lhs, const DMatDMatMapExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( rhs.lhs_ );  // Evaluation of the left-hand side dense matrix operand
      RT B( rhs.rhs_ );  // Evaluation of the right-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == (*lhs).columns()  , "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (*lhs).columns()  , "Invalid number of columns" );

      smpSchurAssign( *lhs, map( A, B, rhs.op_ ) );
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
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER( MT1, MT2 );
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
/*!\brief Elementwise evaluation of the given binary operation on each single element of the
//        dense matrices \a lhs and \a rhs.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix operand.
// \param rhs The right-hand side dense matrix operand.
// \param op The custom, binary operation.
// \return The binary operation applied to each single element of \a lhs and \a rhs.
// \exception std::invalid_argument Matrix sizes do not match.
//
// The \a map() function evaluates the given binary operation on each single element of the input
// matrices \a lhs and \a rhs. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a map() function:

   \code
   blaze::DynamicMatrix<double> A, B, C;
   // ... Resizing and initialization
   C = map( A, B, []( double x, double y ){ return std::min( x, y ); } );
   \endcode

// In case the current number of rows and columns of the two given matrices don't match, a
// \a std::invalid_argument is thrown.
*/
template< typename MT1   // Type of the left-hand side dense matrix
        , typename MT2   // Type of the right-hand side dense matrix
        , bool SO        // Storage order
        , typename OP >  // Type of the custom operation
inline decltype(auto)
   map( const DenseMatrix<MT1,SO>& lhs, const DenseMatrix<MT2,SO>& rhs, OP op )
{
   BLAZE_FUNCTION_TRACE;

   if( (*lhs).rows() != (*rhs).rows() || (*lhs).columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   using ReturnType = const DMatDMatMapExpr<MT1,MT2,OP,SO>;
   return ReturnType( *lhs, *rhs, std::move(op) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Elementwise evaluation of the given ternary operation on each single element of the
//        dense matrices \a dm1, \a dm2 and \a dm3.
// \ingroup dense_matrix
//
// \param dm1 The first dense matrix operand.
// \param dm2 The second dense matrix operand.
// \param dm3 The third dense matrix operand.
// \param op The custom, ternary operation.
// \return The ternary operation applied to each single element of the three matrices.
// \exception std::invalid_argument Matrix sizes do not match.
//
// The \a map() function evaluates the given ternary operation on each single element of the
// input matrices \a dm1, \a dm2, and \a dm3. The function returns an expression representing
// this operation.\n
// In case the current number of rows and columns of the three given matrices don't match, a
// \a std::invalid_argument is thrown.
*/
template< typename MT1   // Type of the first dense matrix
        , typename MT2   // Type of the second dense matrix
        , typename MT3   // Type of the third dense matrix
        , bool SO        // Storage order
        , typename OP >  // Type of the custom operation
inline decltype(auto)
   map( const DenseMatrix<MT1,SO>& dm1, const DenseMatrix<MT2,SO>& dm2,
        const DenseMatrix<MT3,SO>& dm3, OP op )
{
   BLAZE_FUNCTION_TRACE;

   const MakePair mp{};
   return map( map( map( *dm1, *dm2, mp ), *dm3, mp ), join( std::move(op) ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Elementwise evaluation of the given 4-ary operation on each single element of the
//        dense matrices \a dm1, \a dm2, \a dm3 and \a dm4.
// \ingroup dense_matrix
//
// \param dm1 The first dense matrix operand.
// \param dm2 The second dense matrix operand.
// \param dm3 The third dense matrix operand.
// \param dm4 The fourth dense matrix operand.
// \param op The custom, 4-ary operation.
// \return The 4-ary operation applied to each single element of the four matrices.
// \exception std::invalid_argument Matrix sizes do not match.
//
// The \a map() function evaluates the given 4-ary operation on each single element of the input
// matrices \a dm1, \a dm2, \a dm3, and \a dm4. The function returns an expression representing
// this operation.\n
// In case the current number of rows and columns of the four given matrices don't match, a
// \a std::invalid_argument is thrown.
*/
template< typename MT1   // Type of the first dense matrix
        , typename MT2   // Type of the second dense matrix
        , typename MT3   // Type of the third dense matrix
        , typename MT4   // Type of the fourth dense matrix
        , bool SO        // Storage order
        , typename OP >  // Type of the custom operation
inline decltype(auto)
   map( const DenseMatrix<MT1,SO>& dm1, const DenseMatrix<MT2,SO>& dm2,
        const DenseMatrix<MT3,SO>& dm3, const DenseMatrix<MT4,SO>& dm4, OP op )
{
   BLAZE_FUNCTION_TRACE;

   const MakePair mp{};
   return map( map( map( map( *dm1, *dm2, mp ), *dm3, mp ), *dm4, mp ), join( std::move(op) ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Elementwise evaluation of the given 5-ary operation on each single element of the
//        dense matrices \a dm1, \a dm2, \a dm3, \a dm4, and \a dm5.
// \ingroup dense_matrix
//
// \param dm1 The first dense matrix operand.
// \param dm2 The second dense matrix operand.
// \param dm3 The third dense matrix operand.
// \param dm4 The fourth dense matrix operand.
// \param dm5 The fifth dense matrix operand.
// \param op The custom, 5-ary operation.
// \return The 5-ary operation applied to each single element of the five matrices.
// \exception std::invalid_argument Matrix sizes do not match.
//
// The \a map() function evaluates the given 5-ary operation on each single element of the input
// matrices \a dm1, \a dm2, \a dm3, \a dm4, and \a dm5. The function returns an expression
// representing this operation.\n
// In case the current number of rows and columns of the five given matrices don't match, a
// \a std::invalid_argument is thrown.
*/
template< typename MT1   // Type of the first dense matrix
        , typename MT2   // Type of the second dense matrix
        , typename MT3   // Type of the third dense matrix
        , typename MT4   // Type of the fourth dense matrix
        , typename MT5   // Type of the fifth dense matrix
        , bool SO        // Storage order
        , typename OP >  // Type of the custom operation
inline decltype(auto)
   map( const DenseMatrix<MT1,SO>& dm1, const DenseMatrix<MT2,SO>& dm2,
        const DenseMatrix<MT3,SO>& dm3, const DenseMatrix<MT4,SO>& dm4,
        const DenseMatrix<MT5,SO>& dm5, OP op )
{
   BLAZE_FUNCTION_TRACE;

   const MakePair mp{};
   return map( map( map( map( map( *dm1, *dm2, mp ), *dm3, mp ), *dm4, mp ), *dm5, mp ), join( std::move(op) ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Elementwise evaluation of the given 6-ary operation on each single element of the
//        dense matrices \a dm1, \a dm2, \a dm3, \a dm4, \a dm5, and \a dm6.
// \ingroup dense_matrix
//
// \param dm1 The first dense matrix operand.
// \param dm2 The second dense matrix operand.
// \param dm3 The third dense matrix operand.
// \param dm4 The fourth dense matrix operand.
// \param dm5 The fifth dense matrix operand.
// \param dm6 The sixth dense matrix operand.
// \param op The custom, 6-ary operation.
// \return The 6-ary operation applied to each single element of the six matrices.
// \exception std::invalid_argument Matrix sizes do not match.
//
// The \a map() function evaluates the given 6-ary operation on each single element of the input
// matrices \a dm1, \a dm2, \a dm3, \a dm4, \a dm5, and \a dm6. The function returns an expression
// representing this operation.\n
// In case the current number of rows and columns of the six given matrices don't match, a
// \a std::invalid_argument is thrown.
*/
template< typename MT1   // Type of the first dense matrix
        , typename MT2   // Type of the second dense matrix
        , typename MT3   // Type of the third dense matrix
        , typename MT4   // Type of the fourth dense matrix
        , typename MT5   // Type of the fifth dense matrix
        , typename MT6   // Type of the sixth dense matrix
        , bool SO        // Storage order
        , typename OP >  // Type of the custom operation
inline decltype(auto)
   map( const DenseMatrix<MT1,SO>& dm1, const DenseMatrix<MT2,SO>& dm2,
        const DenseMatrix<MT3,SO>& dm3, const DenseMatrix<MT4,SO>& dm4,
        const DenseMatrix<MT5,SO>& dm5, const DenseMatrix<MT6,SO>& dm6, OP op )
{
   BLAZE_FUNCTION_TRACE;

   const MakePair mp{};
   return map( map( map( map( map( map( *dm1, *dm2, mp ), *dm3, mp ), *dm4, mp ), *dm5, mp ), *dm6, mp ), join( std::move(op) ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the componentwise minimum of the dense matrices \a lhs and \a rhs.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix operand.
// \param rhs The right-hand side dense matrix operand.
// \return The resulting dense matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// This function computes the componentwise minimum of the two dense matrices \a lhs and \a rhs.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a min() function:

   \code
   blaze::DynamicMatrix<double> A, B, C;
   // ... Resizing and initialization
   C = min( A, B );
   \endcode

// In case the current number of rows and columns of the two given matrices don't match, a
// \a std::invalid_argument is thrown.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline decltype(auto)
   min( const DenseMatrix<MT1,SO1>& lhs, const DenseMatrix<MT2,SO2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return map( *lhs, *rhs, Min() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the componentwise maximum of the dense matrices \a lhs and \a rhs.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix operand.
// \param rhs The right-hand side dense matrix operand.
// \return The resulting dense matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// This function computes the componentwise maximum of the two dense matrices \a lhs and \a rhs.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a max() function:

   \code
   blaze::DynamicMatrix<double> A, B, C;
   // ... Resizing and initialization
   C = max( A, B );
   \endcode

// In case the current number of rows and columns of the two given matrices don't match, a
// \a std::invalid_argument is thrown.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline decltype(auto)
   max( const DenseMatrix<MT1,SO1>& lhs, const DenseMatrix<MT2,SO2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return map( *lhs, *rhs, Max() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the componentwise hypotenous for the dense matrices \a lhs and \a rhs.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix operand.
// \param rhs The right-hand side dense matrix operand.
// \return The resulting dense matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// The \a hypot() function computes the componentwise hypotenous for the two dense matrices
// \a lhs and \a rhs. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a hypot() function:

   \code
   blaze::DynamicMatrix<double> A, B, C;
   // ... Resizing and initialization
   C = hypot( A, B );
   \endcode

// In case the current number of rows and columns of the two given matrices don't match, a
// \a std::invalid_argument is thrown.

*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline decltype(auto)
   hypot( const DenseMatrix<MT1,SO1>& lhs, const DenseMatrix<MT2,SO2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return map( *lhs, *rhs, Hypot() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the componentwise exponential value for the dense matrices \a lhs and \a rhs.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix operand.
// \param rhs The right-hand side dense matrix operand.
// \return The resulting dense matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// The \a pow() function computes the componentwise exponential value for the two dense matrices
// \a lhs and \a rhs. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a pow() function:

   \code
   blaze::DynamicMatrix<double> A, B, C;
   // ... Resizing and initialization
   C = pow( A, B );
   \endcode

// In case the current number of rows and columns of the two given matrices don't match, a
// \a std::invalid_argument is thrown.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline decltype(auto)
   pow( const DenseMatrix<MT1,SO1>& lhs, const DenseMatrix<MT2,SO2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return map( *lhs, *rhs, Pow() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the multi-valued inverse tangent of the dense matrices \a lhs and \a rhs.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix operand.
// \param rhs The right-hand side dense matrix operand.
// \return The resulting dense matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// This function computes the multi-valued inverse tangent of the two dense matrix \a lhs and
// \a rhs. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a max() function:

   \code
   blaze::DynamicMatrix<double> A, B, C;
   // ... Resizing and initialization
   C = atan2( A, B );
   \endcode

// In case the current number of rows and columns of the two given matrices don't match, a
// \a std::invalid_argument is thrown.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline decltype(auto)
   atan2( const DenseMatrix<MT1,SO1>& lhs, const DenseMatrix<MT2,SO2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return map( *lhs, *rhs, Atan2() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Elementwise conditional selection of values from the dense matrices \a lhs and \a rhs.
// \ingroup dense_matrix
//
// \param cond The dense matrix containing the selection conditions.
// \param lhs The true-case dense matrix.
// \param rhs The false-case dense matrix.
// \return The resulting dense matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// This function performs an elementwise conditional selection of values from the two given dense
// matrices \a lhs and \a rhs. In case an element in the \a cond matrix evaluates to \a true, the
// according element of \a lhs is selected, in case the \a cond element evaluates to \a false, the
// according element of \a rhs is selected. The function returns an expression representing this
// operation.\n
// The following example demonstrates the use of the \a selec() function:

   \code
   blaze::DynamicMatrix<bool> cond{ { true, false }, { true false } };
   blaze::DynamicMatrix<int> A{ { 1, -1 }, { 1, -1 } };
   blaze::DynamicMatrix<int> B{ { -2, 2 }, { -2, 2 } };
   blaze::DynamicMatrix<int> C;
   // ... Resizing and initialization

   C = select( cond, A, B );  // Results in ( 1, 2 ) ( 1, 2 )
   \endcode

// In case the current number of rows and columns of the three given matrices don't match, a
// \a std::invalid_argument is thrown.
*/
template< typename MT1  // Type of the conditional dense matrix
        , typename MT2  // Type of the true-case dense matrix
        , typename MT3  // Type of the false-case dense matrix
        , bool SO >     // Storage order
inline decltype(auto)
   select( const DenseMatrix<MT1,SO>& cond, const DenseMatrix<MT2,SO>& lhs, const DenseMatrix<MT3,SO>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return map( cond, lhs, rhs, []( bool c, const auto& a, const auto& b ) { return c ? a : b; } );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Left-shift operator for the elementwise left-shift of a dense matrix.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix to be shifted.
// \param rhs The right-hand side dense matrix of bits to shift.
// \return The left-shifted dense matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// This operator represents the elementwise left-shift of a given dense matrix:

   \code
   blaze::DynamicMatrix<unsigned int> A, B, C;
   // ... Resizing and initialization
   C = A << B;
   \endcode

// In case the current number of rows and columns of the two given matrices don't match, a
// \a std::invalid_argument is thrown.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline decltype(auto)
   operator<<( const DenseMatrix<MT1,SO1>& lhs, const DenseMatrix<MT2,SO2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return map( *lhs, *rhs, ShiftLV() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Right-shift operator for the elementwise right-shift of a dense matrix.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix to be shifted.
// \param rhs The right-hand side dense matrix of bits to shift.
// \return The right-shifted dense matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// This operator represents the elementwise right-shift of a given dense matrix:

   \code
   blaze::DynamicMatrix<unsigned int> A, B, C;
   // ... Resizing and initialization
   C = A >> B;
   \endcode

// In case the current number of rows and columns of the two given matrices don't match, a
// \a std::invalid_argument is thrown.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline decltype(auto)
   operator>>( const DenseMatrix<MT1,SO1>& lhs, const DenseMatrix<MT2,SO2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return map( *lhs, *rhs, ShiftRV() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise AND operator for two dense matrices.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the bitwise AND operation.
// \param rhs The right-hand side dense matrix for the bitwise AND operation.
// \return The resulting dense matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// This operator represents the bitwise AND of the given two dense matrices:

   \code
   blaze::DynamicMatrix<unsigned int> A, B, C;
   // ... Resizing and initialization
   C = A & B;
   \endcode

// In case the current number of rows and columns of the two given matrices don't match, a
// \a std::invalid_argument is thrown.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline decltype(auto)
   operator&( const DenseMatrix<MT1,SO1>& lhs, const DenseMatrix<MT2,SO2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return map( *lhs, *rhs, Bitand() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise OR operator for two dense matrices.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the bitwise OR operation.
// \param rhs The right-hand side dense matrix for the bitwise OR operation.
// \return The resulting dense matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// This operator represents the bitwise OR of the given two dense matrices:

   \code
   blaze::DynamicMatrix<unsigned int> A, B, C;
   // ... Resizing and initialization
   C = A | B;
   \endcode

// In case the current number of rows and columns of the two given matrices don't match, a
// \a std::invalid_argument is thrown.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline decltype(auto)
   operator|( const DenseMatrix<MT1,SO1>& lhs, const DenseMatrix<MT2,SO2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return map( *lhs, *rhs, Bitor() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise XOR operator for two dense matrices.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the bitwise XOR operation.
// \param rhs The right-hand side dense matrix for the bitwise XOR operation.
// \return The resulting dense matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// This operator represents the bitwise XOR of the given two dense matrices:

   \code
   blaze::DynamicMatrix<unsigned int> A, B, C;
   // ... Resizing and initialization
   C = A ^ B;
   \endcode

// In case the current number of rows and columns of the two given matrices don't match, a
// \a std::invalid_argument is thrown.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline decltype(auto)
   operator^( const DenseMatrix<MT1,SO1>& lhs, const DenseMatrix<MT2,SO2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return map( *lhs, *rhs, Bitxor() );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Left-shift operator for the elementwise left-shift of an elementwise left-shift expression
//        (\f$ A=B\ll C\ll D \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side left-shift expression to be shifted.
// \param rhs The right-hand side dense matrix of bits to shift.
// \return The left-shifted dense matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// This function implements a performance optimized treatment of the elementwise left-shift
// operation on an elementwise left-shift expression.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the middle dense matrix
        , bool SO1      // Storage order of the left-hand side left-shift expression
        , typename MT3  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline decltype(auto)
   operator<<( const DMatDMatMapExpr<MT1,MT2,ShiftLV,SO1>& lhs, const DenseMatrix<MT3,SO2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return map( lhs.leftOperand(), lhs.rightOperand() + (*rhs), lhs.operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Right-shift operator for the elementwise right-shift of an elementwise right-shift
//        expression (\f$ A=B\gg C\gg D \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side right-shift expression to be shifted.
// \param rhs The right-hand side dense matrix of bits to shift.
// \return The right-shifted dense matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// This function implements a performance optimized treatment of the elementwise right-shift
// operation on an elementwise right-shift expression.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the middle dense matrix
        , bool SO1      // Storage order of the left-hand side right-shift expression
        , typename MT3  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline decltype(auto)
   operator>>( const DMatDMatMapExpr<MT1,MT2,ShiftRV,SO1>& lhs, const DenseMatrix<MT3,SO2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return map( lhs.leftOperand(), lhs.rightOperand() + (*rhs), lhs.operation() );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL LOGICAL ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Logical AND operator for two dense matrices.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the logical AND operation.
// \param rhs The right-hand side dense matrix for the logical AND operation.
// \return The resulting dense matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// This operator represents the logical AND of the given two dense matrices:

   \code
   blaze::DynamicMatrix<bool> A, B, C;
   // ... Resizing and initialization
   C = A && B;
   \endcode

// In case the current number of rows and columns of the two given matrices don't match, a
// \a std::invalid_argument is thrown.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline decltype(auto)
   operator&&( const DenseMatrix<MT1,SO1>& lhs, const DenseMatrix<MT2,SO2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return map( *lhs, *rhs, And{} );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Logical OR operator for two dense matrices.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the logical OR operation.
// \param rhs The right-hand side dense matrix for the logical OR operation.
// \return The resulting dense matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// This operator represents the logical OR of the given two dense matrices:

   \code
   blaze::DynamicMatrix<bool> A, B, C;
   // ... Resizing and initialization
   C = A || B;
   \endcode

// In case the current number of rows and columns of the two given matrices don't match, a
// \a std::invalid_argument is thrown.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline decltype(auto)
   operator||( const DenseMatrix<MT1,SO1>& lhs, const DenseMatrix<MT2,SO2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return map( *lhs, *rhs, Or{} );
}
//*************************************************************************************************




//=================================================================================================
//
//  ISALIGNED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, typename OP, bool SO >
struct IsAligned< DMatDMatMapExpr<MT1,MT2,OP,SO> >
   : public BoolConstant< IsAligned_v<MT1> && IsAligned_v<MT2> >
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
template< typename MT1, typename MT2, typename OP, bool SO >
struct IsPadded< DMatDMatMapExpr<MT1,MT2,OP,SO> >
   : public BoolConstant< IsPadded_v<MT1> && IsPadded_v<MT2> && IsPaddingEnabled_v<OP> >
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
