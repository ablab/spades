//=================================================================================================
/*!
//  \file blaze/math/expressions/DMatMapExpr.h
//  \brief Header file for the dense matrix map expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DMATMAPEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DMATMAPEXPR_H_


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
#include <blaze/math/expressions/MatMapExpr.h>
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
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/HasMember.h>
#include <blaze/util/typetraits/IsBuiltin.h>
#include <blaze/util/typetraits/IsNumeric.h>
#include <blaze/util/typetraits/IsSame.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DMATMAPEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for the dense matrix map() function.
// \ingroup dense_matrix_expression
//
// The DMatMapExpr class represents the compile time expression for the evaluation of a custom
// operation on each element of a dense matrix via the map() function.
*/
template< typename MT  // Type of the dense matrix
        , typename OP  // Type of the custom operation
        , bool SO >    // Storage order
class DMatMapExpr
   : public MatMapExpr< DenseMatrix< DMatMapExpr<MT,OP,SO>, SO > >
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   using RT = ResultType_t<MT>;    //!< Result type of the dense matrix expression.
   using OT = OppositeType_t<MT>;  //!< Opposite type of the dense matrix expression.
   using ET = ElementType_t<MT>;   //!< Element type of the dense matrix expression.
   using RN = ReturnType_t<MT>;    //!< Return type of the dense matrix expression.
   //**********************************************************************************************

   //**Serial evaluation strategy******************************************************************
   //! Compilation switch for the serial evaluation strategy of the map expression.
   /*! The \a useAssign compile time constant expression represents a compilation switch for
       the serial evaluation strategy of the map expression. In case the given dense matrix
       expression of type \a MT is a computation expression and requires an intermediate
       evaluation, \a useAssign will be set to 1 and the map expression will be evaluated
       via the \a assign function family. Otherwise \a useAssign will be set to 0 and the
       expression will be evaluated via the subscript operator. */
   static constexpr bool useAssign = ( IsComputation_v<MT> && RequiresEvaluation_v<MT> );

   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT2 >
   static constexpr bool UseAssign_v = useAssign;
   /*! \endcond */
   //**********************************************************************************************

   //**Parallel evaluation strategy****************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! This variable template is a helper for the selection of the parallel evaluation strategy.
       In case either the target matrix or the dense matrix operand is not SMP assignable and
       the matrix operand requires an intermediate evaluation, the variable is set to 1 and the
       expression specific evaluation strategy is selected. Otherwise the variable is set to 0
       and the default strategy is chosen. */
   template< typename MT2 >
   static constexpr bool UseSMPAssign_v =
      ( ( !MT2::smpAssignable || !MT::smpAssignable ) && useAssign );
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this DMatMapExpr instance.
   using This = DMatMapExpr<MT,OP,SO>;

   //! Base type of this DMatMapExpr instance.
   using BaseType = MatMapExpr< DenseMatrix<This,SO> >;

   using ResultType    = MapTrait_t<RT,OP>;            //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<ResultType>;   //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<ResultType>;    //!< Resulting element type.

   //! Return type for expression template evaluations.
   using ReturnType = decltype( std::declval<OP>()( std::declval<RN>() ) );

   //! Data type for composite expression templates.
   using CompositeType = If_t< useAssign, const ResultType, const DMatMapExpr& >;

   //! Composite data type of the dense matrix expression.
   using Operand = If_t< IsExpression_v<MT>, const MT, const MT& >;

   //! Data type of the custom unary operation.
   using Operation = OP;
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

      //! ConstIterator type of the dense matrix expression.
      using IteratorType = ConstIterator_t<MT>;
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor for the ConstIterator class.
      //
      // \param it Iterator to the initial matrix element.
      // \param op The custom unary operation.
      */
      inline BLAZE_DEVICE_CALLABLE ConstIterator( IteratorType it, OP op )
         : it_( it )             // Iterator to the current matrix element
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
      /*!\brief Access to the SIMD elements of the matrix.
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
      IteratorType it_;  //!< Iterator to the current matrix element.
      OP           op_;  //!< The custom unary operation.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled =
      ( MT::simdEnabled &&
        If_t< HasSIMDEnabled_v<OP>, GetSIMDEnabled<OP,ET>, HasLoad<OP> >::value );

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = MT::smpAssignable;
   //**********************************************************************************************

   //**SIMD properties*****************************************************************************
   //! The number of elements packed within a single SIMD element.
   static constexpr size_t SIMDSIZE = SIMDTrait<ElementType>::size;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DMatMapExpr class.
   //
   // \param dm The dense matrix operand of the map expression.
   // \param op The custom unary operation.
   */
   inline DMatMapExpr( const MT& dm, OP op ) noexcept
      : dm_( dm )             // Dense matrix of the map expression
      , op_( std::move(op) )  // The custom unary operation
   {}
   //**********************************************************************************************

   //**Access operator*****************************************************************************
   /*!\brief 2D-access to the matrix elements.
   //
   // \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
   // \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator()( size_t i, size_t j ) const noexcept {
      BLAZE_INTERNAL_ASSERT( i < dm_.rows()   , "Invalid row access index"    );
      BLAZE_INTERNAL_ASSERT( j < dm_.columns(), "Invalid column access index" );
      return op_( dm_(i,j) );
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
      if( i >= dm_.rows() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
      }
      if( j >= dm_.columns() ) {
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
      BLAZE_INTERNAL_ASSERT( i < dm_.rows()   , "Invalid row access index"    );
      BLAZE_INTERNAL_ASSERT( j < dm_.columns(), "Invalid column access index" );
      BLAZE_INTERNAL_ASSERT( !SO || ( i % SIMDSIZE == 0UL ), "Invalid row access index"    );
      BLAZE_INTERNAL_ASSERT( SO  || ( j % SIMDSIZE == 0UL ), "Invalid column access index" );
      return op_.load( dm_.load(i,j) );
   }
   //**********************************************************************************************

   //**Begin function******************************************************************************
   /*!\brief Returns an iterator to the first non-zero element of row/column \a i.
   //
   // \param i The row/column index.
   // \return Iterator to the first non-zero element of row/column \a i.
   */
   inline ConstIterator begin( size_t i ) const {
      return ConstIterator( dm_.begin(i), op_ );
   }
   //**********************************************************************************************

   //**End function********************************************************************************
   /*!\brief Returns an iterator just past the last non-zero element of row/column \a i.
   //
   // \param i The row/column index.
   // \return Iterator just past the last non-zero element of row/column \a i.
   */
   inline ConstIterator end( size_t i ) const {
      return ConstIterator( dm_.end(i), op_ );
   }
   //**********************************************************************************************

   //**Rows function*******************************************************************************
   /*!\brief Returns the current number of rows of the matrix.
   //
   // \return The number of rows of the matrix.
   */
   inline size_t rows() const noexcept {
      return dm_.rows();
   }
   //**********************************************************************************************

   //**Columns function****************************************************************************
   /*!\brief Returns the current number of columns of the matrix.
   //
   // \return The number of columns of the matrix.
   */
   inline size_t columns() const noexcept {
      return dm_.columns();
   }
   //**********************************************************************************************

   //**Operand access******************************************************************************
   /*!\brief Returns the dense matrix operand.
   //
   // \return The dense matrix operand.
   */
   inline Operand operand() const noexcept {
      return dm_;
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
      return IsExpression_v<MT> && dm_.canAlias( alias );
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
      return dm_.isAliased( alias );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the operands of the expression are properly aligned in memory.
   //
   // \return \a true in case the operands are aligned, \a false if not.
   */
   inline bool isAligned() const noexcept {
      return dm_.isAligned();
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const noexcept {
      return dm_.canSMPAssign();
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   Operand   dm_;  //!< Dense matrix of the map expression.
   Operation op_;  //!< The custom unary operation.
   //**********************************************************************************************

   //**Assignment to dense matrices****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense matrix map expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix map
   // expression to a dense matrix. Due to the explicit application of the SFINAE principle,
   // this function can only be selected by the compiler in case the operand requires an
   // intermediate evaluation and the underlying numeric data type of the operand and the
   // target matrix are identical.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order or the target dense matrix
   friend inline auto assign( DenseMatrix<MT2,SO2>& lhs, const DMatMapExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT2> &&
                     IsSame_v< UnderlyingScalar_t<MT>, UnderlyingScalar_t<MT2> > >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      assign( *lhs, rhs.dm_ );
      assign( *lhs, map( *lhs, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to dense matrices****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense matrix map expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix map
   // expression to a dense matrix. Due to the explicit application of the SFINAE principle,
   // this function can only be selected by the compiler in case the operand requires an
   // intermediate evaluation and the underlying numeric data type of the operand and the
   // target vector differ.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order or the target dense matrix
   friend inline auto assign( DenseMatrix<MT2,SO2>& lhs, const DMatMapExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT2> &&
                     !IsSame_v< UnderlyingScalar_t<MT>, UnderlyingScalar_t<MT2> > >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( serial( rhs.dm_ ) );
      assign( *lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse matrices***************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense matrix map expression to a sparse matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side map expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix map
   // expression to a sparse matrix. Due to the explicit application of the SFINAE principle,
   // this operator can only be selected by the compiler in case the operand requires an
   // intermediate evaluation.
   */
   template< typename MT2  // Type of the target sparse matrix
           , bool SO2 >    // Storage order or the target sparse matrix
   friend inline auto assign( SparseMatrix<MT2,SO2>& lhs, const DMatMapExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      using TmpType = If_t< SO == SO2, RT, OT >;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( OT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( OT, !SO );
      BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER( MT2, TmpType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( TmpType );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const TmpType tmp( serial( rhs.dm_ ) );
      assign( *lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense matrices*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a dense matrix map expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense
   // matrix map expression to a dense matrix. Due to the explicit application of the
   // SFINAE principle, this operator can only be selected by the compiler in case the
   // operand requires an intermediate evaluation.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto addAssign( DenseMatrix<MT2,SO2>& lhs, const DMatMapExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( serial( rhs.dm_ ) );
      addAssign( *lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse matrices******************************************************
   // No special implementation for the addition assignment to sparse matrices.
   //**********************************************************************************************

   //**Subtraction assignment to dense matrices****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a dense matrix map expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a dense
   // matrix map expression to a dense matrix. Due to the explicit application of the SFINAE
   // principle, this operator can only be selected by the compiler in case the operand
   // requires an intermediate evaluation.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto subAssign( DenseMatrix<MT2,SO2>& lhs, const DMatMapExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( serial( rhs.dm_ ) );
      subAssign( *lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse matrices***************************************************
   // No special implementation for the subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**Schur product assignment to dense matrices**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Schur product assignment of a dense matrix map expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a dense
   // matrix map expression to a dense matrix. Due to the explicit application of the SFINAE
   // principle, this operator can only be selected by the compiler in case the operand
   // requires an intermediate evaluation.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto schurAssign( DenseMatrix<MT2,SO2>& lhs, const DMatMapExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( serial( rhs.dm_ ) );
      schurAssign( *lhs, map( tmp, rhs.op_ ) );
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
   /*!\brief SMP assignment of a dense matrix map expression to a row-major dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a dense matrix
   // map expression to a row-major dense matrix. Due to the explicit application of the
   // SFINAE principle, this operator can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected and the underlying
   // numeric data type of the operand and the target matrix are identical.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order or the target dense matrix
   friend inline auto smpAssign( DenseMatrix<MT2,SO2>& lhs, const DMatMapExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT2> &&
                     IsSame_v< UnderlyingScalar_t<MT>, UnderlyingScalar_t<MT2> > >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      smpAssign( *lhs, rhs.dm_ );
      smpAssign( *lhs, map( *lhs, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to dense matrices************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a dense matrix map expression to a row-major dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a dense matrix map
   // expression to a row-major dense matrix. Due to the explicit application of the SFINAE
   // principle, this operator can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected and the underlying numeric data type
   // of the operand and the target vector differ.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order or the target dense matrix
   friend inline auto smpAssign( DenseMatrix<MT2,SO2>& lhs, const DMatMapExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT2> &&
                     !IsSame_v< UnderlyingScalar_t<MT>, UnderlyingScalar_t<MT2> > >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( rhs.dm_ );
      smpAssign( *lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to sparse matrices***********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a dense matrix map expression to a sparse matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side map expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a dense matrix map
   // expression to a sparse matrix. Due to the explicit application of the SFINAE principle,
   // this operator can only be selected by the compiler in case the expression specific
   // parallel evaluation strategy is selected.
   */
   template< typename MT2  // Type of the target sparse matrix
           , bool SO2 >    // Storage order or the target sparse matrix
   friend inline auto smpAssign( SparseMatrix<MT2,SO2>& lhs, const DMatMapExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      using TmpType = If_t< SO == SO2, RT, OT >;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( OT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( OT, !SO );
      BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER( MT2, TmpType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( TmpType );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const TmpType tmp( rhs.dm_ );
      smpAssign( *lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to dense matrices***************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a dense matrix map expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a dense
   // matrix map expression to a dense matrix. Due to the explicit application of the SFINAE
   // principle, this operator can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto smpAddAssign( DenseMatrix<MT2,SO2>& lhs, const DMatMapExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( rhs.dm_ );
      smpAddAssign( *lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to sparse matrices**************************************************
   // No special implementation for the SMP addition assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP subtraction assignment to dense matrices************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a dense matrix map expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a
   // dense matrix map expression to a dense matrix. Due to the explicit application of
   // the SFINAE principle, this operator can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto smpSubAssign( DenseMatrix<MT2,SO2>& lhs, const DMatMapExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( rhs.dm_ );
      smpSubAssign( *lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to sparse matrices***********************************************
   // No special implementation for the SMP subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP Schur product assignment to dense matrices**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP Schur product assignment of a dense matrix map expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized SMP Schur product assignment of a
   // dense matrix map expression to a dense matrix. Due to the explicit application of the
   // SFINAE principle, this operator can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto smpSchurAssign( DenseMatrix<MT2,SO2>& lhs, const DMatMapExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( rhs.dm_ );
      smpSchurAssign( *lhs, map( tmp, rhs.op_ ) );
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
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( MT, SO );
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
/*!\brief Evaluates the given custom operation on each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \param op The custom operation.
// \return The custom operation applied to each single element of \a dm.
//
// The \a map() function evaluates the given custom operation on each element of the input
// matrix \a dm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a map() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = map( A, []( double a ){ return std::sqrt( a ); } );
   \endcode
*/
template< typename MT    // Type of the dense matrix
        , bool SO        // Storage order
        , typename OP >  // Type of the custom operation
inline decltype(auto) map( const DenseMatrix<MT,SO>& dm, OP op )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const DMatMapExpr<MT,OP,SO>;
   return ReturnType( *dm, std::move(op) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Evaluates the given custom operation on each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \param op The custom operation.
// \return The custom operation applied to each single element of \a dm.
//
// The \a forEach() function evaluates the given custom operation on each element of the input
// matrix \a dm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a forEach() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = forEach( A, []( double a ){ return std::sqrt( a ); } );
   \endcode
*/
template< typename MT    // Type of the dense matrix
        , bool SO        // Storage order
        , typename OP >  // Type of the custom operation
inline decltype(auto) forEach( const DenseMatrix<MT,SO>& dm, OP op )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, std::move(op) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the componentwise minimum of a dense matrix \a dm and a scalar.
// \ingroup dense_matrix
//
// \param dm The left-hand side dense matrix operand.
// \param scalar The right-hand side scalar value.
// \return The resulting dense matrix.
//
// This operator computes the componentwise minimum of a dense matrix \a dm and a uniform matrix
// represented by the scalar value \a scalar. The function returns an expression representing this
// operation.\n
// The following example demonstrates the use of the \a min() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = min( A, 0.0 );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , typename ST  // Type of the scalar exponent
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
decltype(auto) min( const DenseMatrix<MT,SO>& dm, ST scalar )
{
   BLAZE_FUNCTION_TRACE;

   using ET = ElementType_t<MT>;
   using ScalarType = If_t< IsNumeric_v<ET> && IsNumeric_v<ST>, MapTrait_t<ET,ST,Min>, ST >;
   return map( *dm, bind2nd( Min(), ScalarType( scalar ) ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the componentwise minimum of a scalar and a dense matrix \a dm.
// \ingroup dense_matrix
//
// \param scalar The left-hand side scalar value.
// \param dm The right-hand side dense matrix operand.
// \return The resulting dense matrix.
//
// This operator computes the componentwise minimum of a uniform matrix represented by the scalar
// value \a scalar and a dense matrix \a dm. The function returns an expression representing this
// operation.\n
// The following example demonstrates the use of the \a min() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = min( 0.0, A );
   \endcode
*/
template< typename ST  // Type of the scalar exponent
        , typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
decltype(auto) min( ST scalar, const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   using ET = ElementType_t<MT>;
   using ScalarType = If_t< IsNumeric_v<ST> && IsNumeric_v<ET>, MapTrait_t<ST,ET,Min>, ST >;
   return map( *dm, bind1st( Min(), ScalarType( scalar ) ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the componentwise maximum of a dense matrix \a dm and a scalar.
// \ingroup dense_matrix
//
// \param dm The left-hand side dense matrix operand.
// \param scalar The right-hand side scalar value.
// \return The resulting dense matrix.
//
// This operator computes the componentwise maximum of a dense matrix \a dm and a uniform matrix
// represented by the scalar value \a scalar. The function returns an expression representing this
// operation.\n
// The following example demonstrates the use of the \a max() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = max( A, 0.0 );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , typename ST  // Type of the scalar exponent
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
decltype(auto) max( const DenseMatrix<MT,SO>& dm, ST scalar )
{
   BLAZE_FUNCTION_TRACE;

   using ET = ElementType_t<MT>;
   using ScalarType = If_t< IsNumeric_v<ET> && IsNumeric_v<ST>, MapTrait_t<ET,ST,Max>, ST >;
   return map( *dm, bind2nd( Max(), ScalarType( scalar ) ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the componentwise maximum of a scalar and a dense matrix \a dm.
// \ingroup dense_matrix
//
// \param scalar The left-hand side scalar value.
// \param dm The right-hand side dense matrix operand.
// \return The resulting dense matrix.
//
// This operator computes the componentwise maximum of a uniform matrix represented by the scalar
// value \a scalar and a dense matrix \a dm. The function returns an expression representing this
// operation.\n
// The following example demonstrates the use of the \a max() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = max( 0.0, A );
   \endcode
*/
template< typename ST  // Type of the scalar exponent
        , typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
decltype(auto) max( ST scalar, const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   using ET = ElementType_t<MT>;
   using ScalarType = If_t< IsNumeric_v<ST> && IsNumeric_v<ET>, MapTrait_t<ST,ET,Max>, ST >;
   return map( *dm, bind1st( Max(), ScalarType( scalar ) ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Applies the \a abs() function to each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The resulting dense matrix.
//
// This function applies the \a abs() function to each element of the input matrix \a dm. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a abs() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = abs( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) abs( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, Abs() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Applies the \a sign() function to each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The resulting dense matrix.
//
// This function applies the \a sign() function to each element of the input matrix \a dm. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a sign() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = sign( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) sign( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, Sign() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Applies the \a floor() function to each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The resulting dense matrix.
//
// This function applies the \a floor() function to each element of the input matrix \a dm. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a floor() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = floor( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) floor( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, Floor() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Applies the \a ceil() function to each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The resulting dense matrix.
//
// This function applies the \a ceil() function to each element of the input matrix \a dm. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a ceil() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = ceil( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) ceil( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, Ceil() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Applies the \a trunc() function to each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The resulting dense matrix.
//
// This function applies the \a trunc() function to each element of the input matrix \a dm. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a trunc() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = trunc( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) trunc( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, Trunc() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Applies the \a round() function to each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The resulting dense matrix.
//
// This function applies the \a round() function to each element of the input matrix \a dm. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a round() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = round( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) round( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, Round() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a matrix containing the complex conjugate of each single element of \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The conjugate complex of each single element of \a dm.
//
// The \a conj() function calculates the complex conjugate of each element of the input matrix
// \a dm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a conj() function:

   \code
   blaze::DynamicMatrix< complex<double> > A, B;
   // ... Resizing and initialization
   B = conj( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) conj( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, Conj() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the conjugate transpose matrix of \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The conjugate transpose of \a dm.
//
// The \a ctrans() function returns an expression representing the conjugate transpose (also
// called adjoint matrix, Hermitian conjugate matrix or transjugate matrix) of the given input
// matrix \a dm.\n
// The following example demonstrates the use of the \a ctrans() function:

   \code
   blaze::DynamicMatrix< complex<double> > A, B;
   // ... Resizing and initialization
   B = ctrans( A );
   \endcode

// Note that the \a ctrans() function has the same effect as manually applying the \a conj() and
// \a trans function in any order:

   \code
   B = trans( conj( A ) );  // Computing the conjugate transpose matrix
   B = conj( trans( A ) );  // Computing the conjugate transpose matrix
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) ctrans( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return trans( conj( *dm ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a matrix containing the real part of each single element of \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The real part of each single element of \a dm.
//
// The \a real() function calculates the real part of each element of the input matrix \a dm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a real() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = real( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) real( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, Real() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a matrix containing the imaginary part of each single element of \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The imaginary part of each single element of \a dm.
//
// The \a imag() function calculates the imaginary part of each element of the input matrix \a dm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a imag() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = imag( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) imag( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, Imag() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a matrix containing the phase angle of each single element of \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The phase angle of each single element of \a dm.
//
// The \a arg() function calculates the phase angle of each element of the input matrix \a dm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a arg() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = arg( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) arg( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, Arg() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the square root of each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix; all elements must be in the range \f$[0..\infty)\f$.
// \return The square root of each single element of \a dm.
//
// The \a sqrt() function computes the square root of each element of the input matrix \a dm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a sqrt() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = sqrt( A );
   \endcode

// \note All elements are expected to be in the range \f$[0..\infty)\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) sqrt( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, Sqrt() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the inverse square root of each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix; all elements must be in the range \f$(0..\infty)\f$.
// \return The inverse square root of each single element of \a dm.
//
// The \a invsqrt() function computes the inverse square root of each element of the input matrix
// \a dm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a invsqrt() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = invsqrt( A );
   \endcode

// \note All elements are expected to be in the range \f$(0..\infty)\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) invsqrt( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, InvSqrt() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the cubic root of each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix; all elements must be in the range \f$[0..\infty)\f$.
// \return The cubic root of each single element of \a dm.
//
// The \a cbrt() function computes the cubic root of each element of the input matrix \a dm. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a cbrt() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = cbrt( A );
   \endcode

// \note All elements are expected to be in the range \f$[0..\infty)\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) cbrt( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, Cbrt() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the inverse cubic root of each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix; all elements must be in the range \f$(0..\infty)\f$.
// \return The inverse cubic root of each single element of \a dm.
//
// The \a invcbrt() function computes the inverse cubic root of each element of the input matrix
// \a dm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a invcbrt() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = invcbrt( A );
   \endcode

// \note All elements are expected to be in the range \f$(0..\infty)\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) invcbrt( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, InvCbrt() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Restricts each single element of the dense matrix \a dm to the range \f$[min..max]\f$.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \param min The lower delimiter.
// \param max The upper delimiter.
// \return The matrix with restricted elements.
//
// The \a clamp() function restricts each element of the input matrix \a dm to the range
// \f$[min..max]\f$. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a clamp() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = clamp( A, -1.0, 1.0 );
   \endcode
*/
template< typename MT    // Type of the dense matrix
        , bool SO        // Storage order
        , typename DT >  // Type of the delimiters
inline decltype(auto) clamp( const DenseMatrix<MT,SO>& dm, const DT& min, const DT& max )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, bind2nd( bind3rd( Clamp(), max ), min ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the exponential value for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \param exp The scalar exponent.
// \return The exponential value of each single element of \a dm.
//
// The \a pow() function computes the exponential value for each element of the input matrix
// \a dm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a pow() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = pow( A, 4.2 );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO      // Storage order
        , typename ST  // Type of the scalar exponent
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
inline decltype(auto) pow( const DenseMatrix<MT,SO>& dm, ST exp )
{
   BLAZE_FUNCTION_TRACE;

   using ET = ElementType_t<MT>;
   using ScalarType = If_t< IsNumeric_v<ET> && IsNumeric_v<ST>, MultTrait_t<ET,ST>, ST >;
   return map( *dm, blaze::bind2nd( Pow(), ScalarType( exp ) ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes \f$ e^x \f$ for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The resulting dense matrix.
//
// The \a exp() function computes \f$ e^x \f$ for each element of the input matrix \a dm. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a exp() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = exp( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) exp( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, Exp() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes \f$ 2^x \f$ for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The resulting dense matrix.
//
// The \a exp2() function computes \f$ 2^x \f$ for each element of the input matrix \a dm. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a exp2() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = exp2( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) exp2( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, Exp2() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes \f$ 10^x \f$ for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The resulting dense matrix.
//
// The \a exp10() function computes \f$ 10^x \f$ for each element of the input matrix \a dm. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a exp10() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = exp10( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) exp10( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, Exp10() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the natural logarithm for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix; all elements must be in the range \f$[0..\infty)\f$.
// \return The natural logarithm of each single element of \a dm.
//
// The \a log() function computes natural logarithm for each element of the input matrix \a dm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a log() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = log( A );
   \endcode

// \note All elements are expected to be in the range \f$[0..\infty)\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) log( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, Log() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the binary logarithm for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix; all elements must be in the range \f$[0..\infty)\f$.
// \return The binary logarithm of each single element of \a dm.
//
// The \a log2() function computes binary logarithm for each element of the input matrix \a dm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a log2() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = log2( A );
   \endcode

// \note All elements are expected to be in the range \f$[0..\infty)\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) log2( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, Log2() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the common logarithm for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix; all elements must be in the range \f$[0..\infty)\f$.
// \return The common logarithm of each single element of \a dm.
//
// The \a log10() function computes common logarithm for each element of the input matrix \a dm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a log10() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = log10( A );
   \endcode

// \note All elements are expected to be in the range \f$[0..\infty)\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) log10( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, Log10() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the natural logarithm of x+1 for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix; all elements must be in the range \f$[0..\infty)\f$.
// \return The natural logarithm of x+1 for each single element of \a dm.
//
// The \a log1p() function computes the natural logarithm of x+1 for each element of the input
// matrix \a dm. This may be preferred over the natural logarithm for higher precision computing
// the natural logarithm of a quantity very close to 1. The function returns an expression
// representing this operation.\n
// The following example demonstrates the use of the \a log1p() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = log1p( A );
   \endcode

// \note All elements are expected to be in the range \f$(-1..\infty)\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) log1p( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, Log1p() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the natural logarithm of the absolute value of the gamma function for each
//        single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix; all elements must be in the range \f$[0..\infty)\f$.
// \return The natural logarithm of the absolute value of the gamma function of each single element of \a dm.
//
// The \a lgamma() function computes the natural logarithm of the absolute value of the gamma
// function for each element of the input matrix \a dm. The function returns an expression
// representing this operation.\n
// The following example demonstrates the use of the \a lgamma() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = lgamma( A );
   \endcode

// \note All elements are expected to be in the range \f$[0..\infty)\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) lgamma( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, LGamma() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the sine for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The sine of each single element of \a dm.
//
// The \a sin() function computes the sine for each element of the input matrix \a dm. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a sin() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = sin( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) sin( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, Sin() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the inverse sine for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix; all elements must be in the range \f$[-1..1]\f$.
// \return The inverse sine of each single element of \a dm.
//
// The \a asin() function computes the inverse sine for each element of the input matrix \a dm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a asin() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = asin( A );
   \endcode

// \note All elements are expected to be in the range \f$[-1..1]\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) asin( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, Asin() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the hyperbolic sine for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The hyperbolic sine of each single element of \a dm.
//
// The \a sinh() function computes the hyperbolic sine for each element of the input matrix \a dm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a sinh() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = sinh( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) sinh( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, Sinh() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the inverse hyperbolic sine for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The inverse hyperbolic sine of each single element of \a dm.
//
// The \a asinh() function computes the inverse hyperbolic sine for each element of the input
// matrix \a dm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a asinh() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = asinh( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) asinh( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, Asinh() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the cosine for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The cosine of each single element of \a dm.
//
// The \a cos() function computes the cosine for each element of the input matrix \a dm. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a cos() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = cos( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) cos( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, Cos() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the inverse cosine for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix; all elements must be in the range \f$[-1..1]\f$.
// \return The inverse cosine of each single element of \a dm.
//
// The \a acos() function computes the inverse cosine for each element of the input matrix \a dm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a acos() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = acos( A );
   \endcode

// \note All elements are expected to be in the range \f$[-1..1]\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) acos( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, Acos() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the hyperbolic cosine for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The hyperbolic cosine of each single element of \a dm.
//
// The \a cosh() function computes the hyperbolic cosine for each element of the input matrix
// \a dm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a cosh() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = cosh( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) cosh( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, Cosh() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the inverse hyperbolic cosine for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix; all elements must be in the range \f$[1..\infty)\f$.
// \return The inverse hyperbolic cosine of each single element of \a dm.
//
// The \a acosh() function computes the inverse hyperbolic cosine for each element of the input
// matrix \a dm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a acosh() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = acosh( A );
   \endcode

// \note All elements are expected to be in the range \f$[1..\infty)\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) acosh( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, Acosh() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the tangent for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The tangent of each single element of \a dm.
//
// The \a tan() function computes the tangent for each element of the input matrix \a dm. The
// function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a tan() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = tan( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) tan( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, Tan() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the inverse tangent for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The inverse tangent of each single element of \a dm.
//
// The \a atan() function computes the inverse tangent for each element of the input matrix \a dm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a atan() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = atan( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) atan( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, Atan() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the hyperbolic tangent for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix; all elements must be in the range \f$[-1..1]\f$.
// \return The hyperbolic tangent of each single element of \a dm.
//
// The \a tanh() function computes the hyperbolic tangent for each element of the input matrix
// \a dm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a tanh() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = tanh( A );
   \endcode

// \note All elements are expected to be in the range \f$[-1..1]\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) tanh( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, Tanh() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the inverse hyperbolic tangent for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix; all elements must be in the range \f$[-1..1]\f$.
// \return The inverse hyperbolic tangent of each single element of \a dm.
//
// The \a atanh() function computes the inverse hyperbolic tangent for each element of the input
// matrix \a dm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a atanh() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = atanh( A );
   \endcode

// \note All elements are expected to be in the range \f$[-1..1]\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) atanh( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, Atanh() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the error function for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The error function of each single element of \a dm.
//
// The \a erf() function computes the error function for each element of the input matrix \a dm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a erf() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = erf( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) erf( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, Erf() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the complementary error function for each single element of the dense matrix \a dm.
// \ingroup dense_matrix
//
// \param dm The input matrix.
// \return The complementary error function of each single element of \a dm.
//
// The \a erfc() function computes the complementary error function for each element of the input
// matrix \a dm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a erfc() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = erfc( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) erfc( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *dm, Erfc() );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Absolute value function for dense matrix absolute value expressions.
// \ingroup dense_matrix
//
// \param dm The absolute value dense matrix expression.
// \return The absolute value of each single element of \a dm.
//
// This function implements a performance optimized treatment of the absolute value operation
// on a dense matrix absolute value expression.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) abs( const DMatMapExpr<MT,Abs,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return dm;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Applies the \a sign() function for dense matrix \a sign() expressions.
// \ingroup dense_matrix
//
// \param dm The dense matrix \a sign() expression.
// \return The sign of each single element of \a dm.
//
// This function implements a performance optimized treatment of the \a sign() operation on a
// dense matrix \a sign() expression.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) sign( const DMatMapExpr<MT,Sign,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return dm;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Applies the \a floor() function to a dense matrix \a floor() expressions.
// \ingroup dense_matrix
//
// \param dm The dense matrix \a floor() expression.
// \return The resulting dense matrix.
//
// This function implements a performance optimized treatment of the \a floor() operation on
// a dense matrix \a floor() expression.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) floor( const DMatMapExpr<MT,Floor,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return dm;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Applies the \a ceil() function to a dense matrix \a ceil() expressions.
// \ingroup dense_matrix
//
// \param dm The dense matrix \a ceil() expression.
// \return The resulting dense matrix.
//
// This function implements a performance optimized treatment of the \a ceil() operation on
// a dense matrix \a ceil() expression.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) ceil( const DMatMapExpr<MT,Ceil,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return dm;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Applies the \a trunc() function to a dense matrix \a trunc() expressions.
// \ingroup dense_matrix
//
// \param dm The dense matrix \a trunc() expression.
// \return The resulting dense matrix.
//
// This function implements a performance optimized treatment of the \a trunc() operation on
// a dense matrix \a trunc() expression.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) trunc( const DMatMapExpr<MT,Trunc,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return dm;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Applies the \a round() function to a dense matrix \a round() expressions.
// \ingroup dense_matrix
//
// \param dm The dense matrix \a round() expression.
// \return The resulting dense matrix.
//
// This function implements a performance optimized treatment of the \a round() operation on
// a dense matrix \a round() expression.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) round( const DMatMapExpr<MT,Round,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return dm;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Complex conjugate function for complex conjugate dense matrix expressions.
// \ingroup dense_matrix
//
// \param dm The complex conjugate dense matrix expression.
// \return The original dense matrix.
//
// This function implements a performance optimized treatment of the complex conjugate operation
// on a dense matrix complex conjugate expression. It returns an expression representing the
// original dense matrix:

   \code
   blaze::DynamicMatrix< complex<double> > A, B;
   // ... Resizing and initialization
   B = conj( conj( A ) );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) conj( const DMatMapExpr<MT,Conj,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return dm.operand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Complex conjugate function for conjugate transpose dense matrix expressions.
// \ingroup dense_matrix
//
// \param dm The conjugate transpose dense matrix expression.
// \return The transpose dense matrix.
//
// This function implements a performance optimized treatment of the complex conjugate operation
// on a dense matrix conjugate transpose expression. It returns an expression representing the
// transpose of the dense matrix:

   \code
   blaze::DynamicMatrix< complex<double> > A, B;
   // ... Resizing and initialization
   B = conj( ctrans( A ) );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) conj( const DMatTransExpr<DMatMapExpr<MT,Conj,SO>,!SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return trans( dm.operand().operand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Real part function for real part dense matrix expressions.
// \ingroup dense_matrix
//
// \param dm The real part dense matrix expression.
// \return The real part of each single element of \a dm.
//
// This function implements a performance optimized treatment of the real part operation on
// a dense matrix real part expression.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) real( const DMatMapExpr<MT,Real,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return dm;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Addition operator for the addition of a dense matrix and a scalar value (\f$ A=B+s \f$).
// \ingroup dense_matrix
//
// \param mat The left-hand side dense matrix for the addition.
// \param scalar The right-hand side scalar value for the addition.
// \return The matrix sum.
//
// This operator represents the elementwise addition of a dense matrix and a uniform matrix
// represented by a scalar value:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = A + 1.25;
   \endcode

// The operator returns an expression representing a dense matrix of the higher-order element
// type of the involved data types \a MT::ElementType and \a ST. Note that this operator only
// works for scalar values of built-in data type.
*/
template< typename MT  // Type of the left-hand side dense matrix
        , bool SO      // Storage order of the left-hand side dense matrix
        , typename ST  // Type of the right-hand side scalar
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
inline decltype(auto) operator+( const DenseMatrix<MT,SO>& mat, ST scalar )
{
   BLAZE_FUNCTION_TRACE;

   using ET = ElementType_t<MT>;
   using ScalarType = If_t< IsNumeric_v<ET> && IsNumeric_v<ST>, AddTrait_t<ET,ST>, ST >;
   return map( *mat, blaze::bind2nd( Add{}, ScalarType( scalar ) ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition operator for the addition of a scalar value and a dense matrix (\f$ A=s+B \f$).
// \ingroup dense_matrix
//
// \param scalar The left-hand side scalar value for the addition.
// \param mat The right-hand side dense matrix for the addition.
// \return The matrix sum.
//
// This operator represents the elementwise addition of a uniform matrix represented by a scalar
// value and a dense matrix:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = 1.25 + A;
   \endcode

// The operator returns an expression representing a dense matrix of the higher-order element
// type of the involved data types \a MT::ElementType and \a ST. Note that this operator only
// works for scalar values of built-in data type.
*/
template< typename ST  // Type of the left-hand side scalar
        , typename MT  // Type of the right-hand side dense matrix
        , bool SO      // Storage order of the right-hand side dense matrix
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
inline decltype(auto) operator+( ST scalar, const DenseMatrix<MT,SO>& mat )
{
   BLAZE_FUNCTION_TRACE;

   using ET = ElementType_t<MT>;
   using ScalarType = If_t< IsNumeric_v<ST> && IsNumeric_v<ET>, AddTrait_t<ST,ET>, ST >;
   return map( *mat, blaze::bind1st( Add{}, ScalarType( scalar ) ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction operator for the subtraction of a dense matrix and a scalar value
//        (\f$ A=B-s \f$).
// \ingroup dense_matrix
//
// \param mat The left-hand side dense matrix for the subtraction.
// \param scalar The right-hand side scalar value for the subtraction.
// \return The matrix difference.
//
// This operator represents the elementwise subtraction of a uniform matrix represented by a
// scalar value from a dense matrix:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = A - 1.25;
   \endcode

// The operator returns an expression representing a dense matrix of the higher-order element
// type of the involved data types \a MT::ElementType and \a ST. Note that this operator only
// works for scalar values of built-in data type.
*/
template< typename MT  // Type of the left-hand side dense matrix
        , bool SO      // Storage order of the left-hand side dense matrix
        , typename ST  // Type of the right-hand side scalar
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
inline decltype(auto) operator-( const DenseMatrix<MT,SO>& mat, ST scalar )
{
   BLAZE_FUNCTION_TRACE;

   using ET = ElementType_t<MT>;
   using ScalarType = If_t< IsNumeric_v<ET> && IsNumeric_v<ST>, SubTrait_t<ET,ST>, ST >;
   return map( *mat, blaze::bind2nd( Sub{}, ScalarType( scalar ) ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction operator for the subtraction of a scalar value and a dense matrix
//        (\f$ A=s-B \f$).
// \ingroup dense_matrix
//
// \param scalar The left-hand side scalar value for the subtraction.
// \param mat The right-hand side dense matrix for the subtraction.
// \return The matrix difference.
//
// This operator represents the elementwise subtraction of a dense matrix from a uniform matrix
// represented by a scalar value:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = 1.25 - A;
   \endcode

// The operator returns an expression representing a dense matrix of the higher-order element
// type of the involved data types \a MT::ElementType and \a ST. Note that this operator only
// works for scalar values of built-in data type.
*/
template< typename ST  // Type of the left-hand side scalar
        , typename MT  // Type of the right-hand side dense matrix
        , bool SO      // Storage order of the right-hand side dense matrix
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
inline decltype(auto) operator-( ST scalar, const DenseMatrix<MT,SO>& mat )
{
   BLAZE_FUNCTION_TRACE;

   using ET = ElementType_t<MT>;
   using ScalarType = If_t< IsNumeric_v<ST> && IsNumeric_v<ET>, SubTrait_t<ST,ET>, ST >;
   return map( *mat, blaze::bind1st( Sub{}, ScalarType( scalar ) ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division operator for the division of a scalar value and a dense matrix (\f$ A=s/B \f$).
// \ingroup dense_matrix
//
// \param scalar The left-hand side scalar value for the division.
// \param mat The right-hand side dense matrix for the division.
// \return The matrix quotient.
//
// This operator represents the elementwise division of a uniform matrix represented by a scalar
// value and a dense matrix:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = 1.25 / A;
   \endcode

// The operator returns an expression representing a dense matrix of the higher-order element
// type of the involved data types \a MT::ElementType and \a ST. Note that this operator only
// works for scalar values of built-in data type.
*/
template< typename ST  // Type of the left-hand side scalar
        , typename MT  // Type of the right-hand side dense matrix
        , bool SO      // Storage order of the right-hand side dense matrix
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
inline decltype(auto) operator/( ST scalar, const DenseMatrix<MT,SO>& mat )
{
   BLAZE_FUNCTION_TRACE;

   using ET = ElementType_t<MT>;
   using ScalarType = If_t< IsNumeric_v<ST> && IsNumeric_v<ET>, DivTrait_t<ST,ET>, ST >;
   return map( *mat, blaze::bind1st( Div{}, ScalarType( scalar ) ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Left-shift operator for the uniform left-shift of a dense matrix.
// \ingroup dense_matrix
//
// \param mat The dense matrix for the uniform left-shift operation.
// \param count The number of bits to shift all matrix elements.
// \return The resulting matrix.
//
// This operator represents the uniform left-shift of all elements of a dense matrix:

   \code
   blaze::DynamicMatrix<unsigned int> A, B;
   // ... Resizing and initialization
   B = A << 3;
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Transpose flag
inline decltype(auto) operator<<( const DenseMatrix<MT,SO>& mat, int count )
{
   BLAZE_FUNCTION_TRACE;

   return map( *mat, ShiftLI( count ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Right-shift operator for the uniform right-shift of a dense matrix.
// \ingroup dense_matrix
//
// \param mat The dense matrix for the uniform right-shift operation.
// \param count The number of bits to shift all matrix elements.
// \return The resulting matrix.
//
// This operator represents the uniform right-shift of all elements of a dense matrix:

   \code
   blaze::DynamicMatrix<unsigned int> A, B;
   // ... Resizing and initialization
   B = A >> 3;
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Transpose flag
inline decltype(auto) operator>>( const DenseMatrix<MT,SO>& mat, int count )
{
   BLAZE_FUNCTION_TRACE;

   return map( *mat, ShiftRI( count ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise AND operator for the bitwise AND of a dense matrix and a scalar value.
// \ingroup dense_matrix
//
// \param mat The left-hand side dense matrix for the bitwise AND.
// \param scalar The right-hand side scalar value for the bitwise AND.
// \return The resulting matrix.
//
// This operator represents the bitwise AND of a scalar value with all elements of a dense matrix:

   \code
   blaze::DynamicMatrix<unsigned int> A, B;
   // ... Resizing and initialization
   B = A & 7U;
   \endcode
*/
template< typename MT  // Type of the left-hand side dense matrix
        , bool SO      // Storage order of the left-hand side dense matrix
        , typename ST  // Type of the right-hand side scalar
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
inline decltype(auto) operator&( const DenseMatrix<MT,SO>& mat, ST scalar )
{
   BLAZE_FUNCTION_TRACE;

   return map( *mat, blaze::bind2nd( Bitand{}, scalar ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise OR operator for the bitwise OR of a dense matrix and a scalar value.
// \ingroup dense_matrix
//
// \param mat The left-hand side dense matrix for the bitwise OR.
// \param scalar The right-hand side scalar value for the bitwise OR.
// \return The resulting matrix.
//
// This operator represents the bitwise OR of a scalar value with all elements of a dense matrix:

   \code
   blaze::DynamicMatrix<unsigned int> A, B;
   // ... Resizing and initialization
   B = A | 7U;
   \endcode
*/
template< typename MT  // Type of the left-hand side dense matrix
        , bool SO      // Storage order of the left-hand side dense matrix
        , typename ST  // Type of the right-hand side scalar
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
inline decltype(auto) operator|( const DenseMatrix<MT,SO>& mat, ST scalar )
{
   BLAZE_FUNCTION_TRACE;

   return map( *mat, blaze::bind2nd( Bitor{}, scalar ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise XOR operator for the bitwise XOR of a dense matrix and a scalar value.
// \ingroup dense_matrix
//
// \param mat The left-hand side dense matrix for the bitwise XOR.
// \param scalar The right-hand side scalar value for the bitwise XOR.
// \return The resulting matrix.
//
// This operator represents the bitwise XOR of a scalar value with all elements of a dense matrix:

   \code
   blaze::DynamicMatrix<unsigned int> A, B;
   // ... Resizing and initialization
   B = A ^ 7U;
   \endcode
*/
template< typename MT  // Type of the left-hand side dense matrix
        , bool SO      // Storage order of the left-hand side dense matrix
        , typename ST  // Type of the right-hand side scalar
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
inline decltype(auto) operator^( const DenseMatrix<MT,SO>& mat, ST scalar )
{
   BLAZE_FUNCTION_TRACE;

   return map( *mat, blaze::bind2nd( Bitxor{}, scalar ) );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL LOGICAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Logical NOT operator for the logical NOT of a dense matrix.
// \ingroup dense_matrix
//
// \param mat The dense matrix for the logical NOT.
// \return The negated matrix.
//
// This operator represents the logical NOT of all elements of a dense matrix:

   \code
   blaze::DynamicMatrix<bool> A, B;
   // ... Resizing and initialization
   B = !A;
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline decltype(auto) operator!( const DenseMatrix<MT,SO>& mat )
{
   BLAZE_FUNCTION_TRACE;

   return map( *mat, Not{} );
}
//*************************************************************************************************




//=================================================================================================
//
//  ISALIGNED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, typename OP, bool SO >
struct IsAligned< DMatMapExpr<MT,OP,SO> >
   : public IsAligned<MT>
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
template< typename MT, typename OP, bool SO >
struct IsPadded< DMatMapExpr<MT,OP,SO> >
   : public BoolConstant< IsPadded_v<MT> && IsPaddingEnabled_v<OP> >
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
