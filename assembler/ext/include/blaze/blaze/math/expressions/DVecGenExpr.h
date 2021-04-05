//=================================================================================================
/*!
//  \file blaze/math/expressions/DVecGenExpr.h
//  \brief Header file for the dense vector generator expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DVECGENEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DVECGENEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <utility>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/Expression.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/dense/Forward.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/VecGenExpr.h>
#include <blaze/math/shims/Evaluate.h>
#include <blaze/math/shims/Exp10.h>
#include <blaze/math/typetraits/Size.h>
#include <blaze/math/typetraits/TransposeFlag.h>
#include <blaze/math/typetraits/UnderlyingBuiltin.h>
#include <blaze/math/views/Check.h>
#include <blaze/system/HostDevice.h>
#include <blaze/system/TransposeFlag.h>
#include <blaze/util/Assert.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/SmallArray.h>
#include <blaze/util/StaticAssert.h>
#include <blaze/util/TypeList.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsFloatingPoint.h>
#include <blaze/util/typetraits/RemoveCVRef.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DVECGENEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for the dense vector generate() function.
// \ingroup dense_vector_expression
//
// The DVecGenExpr class represents the compile time expression for the generation of a custom
// dense vector via the generate() function.
*/
template< typename VT  // Type of the dense vector
        , typename OP  // Type of the custom unary operation
        , bool TF >    // Transpose flag
class DVecGenExpr
   : public VecGenExpr< DenseVector< DVecGenExpr<VT,OP,TF>, TF > >
   , private Computation
{
 public:
   //**Type definitions****************************************************************************
   //! Type of this DVecGenExpr instance.
   using This = DVecGenExpr<VT,OP,TF>;

   //! Base type of this DVecGenExpr instance.
   using BaseType = VecGenExpr< DenseVector<This,TF> >;

   using ResultType    = RemoveCVRef_t<VT>;    //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<VT>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<VT>;    //!< Resulting element type.

   //! Return type for expression template evaluations.
   using ReturnType = decltype( std::declval<OP>()( std::declval<size_t>() ) );

   //! Data type for composite expression templates.
   using CompositeType = const DVecGenExpr&;

   //! Data type of the custom unary operation.
   using Operation = OP;
   //**********************************************************************************************

   //**ConstIterator class definition**************************************************************
   /*!\brief Iterator over the elements of the dense vector generator expression.
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
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor for the ConstIterator class.
      //
      // \param index Index to the initial vector element.
      // \param op The custom unary operation.
      */
      inline BLAZE_DEVICE_CALLABLE ConstIterator( size_t index, OP op )
         : index_( index )          // Index of the current vector element
         , op_   ( std::move(op) )  // The custom unary operation
      {}
      //*******************************************************************************************

      //**Addition assignment operator*************************************************************
      /*!\brief Addition assignment operator.
      //
      // \param inc The increment of the iterator.
      // \return The incremented iterator.
      */
      inline BLAZE_DEVICE_CALLABLE ConstIterator& operator+=( size_t inc ) {
         index_ += inc;
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
         index_ -= dec;
         return *this;
      }
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented iterator.
      */
      inline BLAZE_DEVICE_CALLABLE ConstIterator& operator++() {
         ++index_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix increment operator***************************************************************
      /*!\brief Post-increment operator.
      //
      // \return The previous position of the iterator.
      */
      inline BLAZE_DEVICE_CALLABLE const ConstIterator operator++( int ) {
         return ConstIterator( index_++, op_ );
      }
      //*******************************************************************************************

      //**Prefix decrement operator****************************************************************
      /*!\brief Pre-decrement operator.
      //
      // \return Reference to the decremented iterator.
      */
      inline BLAZE_DEVICE_CALLABLE ConstIterator& operator--() {
         --index_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix decrement operator***************************************************************
      /*!\brief Post-decrement operator.
      //
      // \return The previous position of the iterator.
      */
      inline BLAZE_DEVICE_CALLABLE const ConstIterator operator--( int ) {
         return ConstIterator( index_--, op_ );
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the element at the current iterator position.
      //
      // \return The resulting value.
      */
      inline BLAZE_DEVICE_CALLABLE ReturnType operator*() const {
         return op_( index_ );
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator==( const ConstIterator& rhs ) const noexcept {
         return index_ == rhs.index_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator!=( const ConstIterator& rhs ) const noexcept {
         return index_ != rhs.index_;
      }
      //*******************************************************************************************

      //**Less-than operator***********************************************************************
      /*!\brief Less-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller, \a false if not.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator<( const ConstIterator& rhs ) const noexcept {
         return index_ < rhs.index_;
      }
      //*******************************************************************************************

      //**Greater-than operator********************************************************************
      /*!\brief Greater-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater, \a false if not.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator>( const ConstIterator& rhs ) const noexcept {
         return index_ > rhs.index_;
      }
      //*******************************************************************************************

      //**Less-or-equal-than operator**************************************************************
      /*!\brief Less-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller or equal, \a false if not.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator<=( const ConstIterator& rhs ) const noexcept {
         return index_ <= rhs.index_;
      }
      //*******************************************************************************************

      //**Greater-or-equal-than operator***********************************************************
      /*!\brief Greater-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater or equal, \a false if not.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator>=( const ConstIterator& rhs ) const noexcept {
         return index_ >= rhs.index_;
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two iterators.
      //
      // \param rhs The right-hand side iterator.
      // \return The number of elements between the two iterators.
      */
      inline BLAZE_DEVICE_CALLABLE DifferenceType operator-( const ConstIterator& rhs ) const noexcept {
         return index_ - rhs.index_;
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
         return ConstIterator( it.index_ + inc, it.op_ );
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
         return ConstIterator( it.index_ + inc, it.op_ );
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
         return ConstIterator( it.index_ - dec, it.op_ );
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      size_t index_;  //!< Index of the current vector element.
      OP     op_;     //!< The custom unary operation.
      //*******************************************************************************************
   };
   //**********************************************************************************************

 public:
   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = false;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = true;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DVecGenExpr class.
   //
   // \param size The size/dimension of the dense vector generator.
   // \param op The custom unary operation.
   */
   inline DVecGenExpr( size_t size, OP&& op ) noexcept
      : size_( size )             // The size/dimension of the dense vector generator
      , op_  ( std::move( op ) )  // The custom unary operation
   {}
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      return op_( index );
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
      if( index >= size_ ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid vector access index" );
      }
      return (*this)[index];
   }
   //**********************************************************************************************

   //**Begin function******************************************************************************
   /*!\brief Returns an iterator to the first non-zero element of the dense vector.
   //
   // \return Iterator to the first non-zero element of the dense vector.
   */
   inline ConstIterator begin() const {
      return ConstIterator( 0UL, op_ );
   }
   //**********************************************************************************************

   //**End function********************************************************************************
   /*!\brief Returns an iterator just past the last non-zero element of the dense vector.
   //
   // \return Iterator just past the last non-zero element of the dense vector.
   */
   inline ConstIterator end() const {
      return ConstIterator( size_, op_ );
   }
   //**********************************************************************************************

   //**Size function*******************************************************************************
   /*!\brief Returns the current size/dimension of the vector.
   //
   // \return The size of the vector.
   */
   inline size_t size() const noexcept {
      return size_;
   }
   //**********************************************************************************************

   //**Operation access****************************************************************************
   /*!\brief Returns a copy of the custom unary operation.
   //
   // \return A copy of the custom unary operation.
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
      MAYBE_UNUSED( alias );
      return false;
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
      MAYBE_UNUSED( alias );
      return false;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the operands of the expression are properly aligned in memory.
   //
   // \return \a true in case the operands are aligned, \a false if not.
   */
   inline bool isAligned() const noexcept {
      return true;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const noexcept {
      return true;
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   size_t size_;   //!< The size/dimension of the dense vector generator.
   Operation op_;  //!< The custom unary operation.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT, TF );
   BLAZE_CONSTRAINT_MUST_NOT_BE_EXPRESSION_TYPE( VT );
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
/*! \cond BLAZE_INTERNAL */
/*!\brief Generates a new dense vector filled via the given custom unary operation.
// \ingroup dense_vector
//
// \param size The size/dimension of the vector.
// \param op The custom unary operation.
// \return The newly generated dense vector.
// \exception std::invalid_argument Invalid size specification.
//
// The \a generate() function returns a dense vector filled elementwise via the given custom unary
// operation. The function returns an expression representing an accordingly filled vector of
// the specified type \a VT.\n
// The following example demonstrates the use of the \a generate() function:

   \code
   using blaze::columnVector;

   // Generates the uniform integer vector ( 2, 2, 2, 2, 2 )
   using A = blaze::DynamicVector<int,columnVector>;
   const auto a = generate<A>( 5UL, []( size_t index ){ return 2; } );

   // Generates the linearly spaced float vector ( 2.1, 3.2, 4.3, 5.4 )
   using B = blaze::DynamicVector<float,columnVector>;
   const auto b = generate<B>( 4UL, []( size_t index ){ return 2.1F + 1.1F*index; } );

   // Generates the logarithmically spaced double vector ( 1.0, 10.0, 100.0, 1000.0 )
   using C = blaze::DynamicVector<double,columnVector>;
   const auto c = generate<C>( 4UL, []( size_t index ){ return blaze::exp10( 1.0 + 1.0*index ); } );

   // Generates the vector of integer vectors ( ( 1, 2 ), ( 2, 3 ), ( 3, 4 ), ( 4, 5 ) )
   using VT = StaticVector<int,2UL>;
   using D = blaze::DynamicVector<VT,columnVector>;
   const auto d = generate<D>( 4UL, []( size_t index ) { return evaluate( VT{ 1, 2 } + index ); } );
   \endcode

// In case the specified size does not match the size of the given vector type \a VT, a
// \a std::invalid_argument is thrown.
*/
template< typename VT    // Type of the dense vector
        , typename OP >  // Type of the custom unary operation
inline decltype(auto) generate( size_t size, OP op )
{
   BLAZE_FUNCTION_TRACE;

   if( Size_v<VT,0UL> != DefaultSize_v && size_t( Size_v<VT,0UL> ) != size ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid size specification" );
   }

   using ReturnType = const DVecGenExpr< VT, OP, TransposeFlag_v<VT> >;
   return ReturnType( size, std::move( op ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Generates a new dense vector filled via the given custom unary operation.
// \ingroup dense_vector
//
// \param size The size/dimension of the vector.
// \param op The custom unary operation.
// \return The newly generated dense vector.
//
// The \c generate() function returns a dense vector filled elementwise via the given custom unary
// operation. By default, the returned vector is a column vector, but this setting can be changed
// via the \c BLAZE_DEFAULT_TRANSPOSE_FLAG switch. Alternatively it is possible to specify the
// transpose flag explicitly.\n
// The following example demonstrates the use of the \a generate() function:

   \code
   using blaze::generate;
   using blaze::columnVector;
   using blaze::rowVector>

   // Generates the uniform integer vector ( 2, 2, 2, 2, 2 )
   blaze::DynamicVector<int,columnVector> a;
   a = generate( 5UL, []( size_t index ){ return 2; } );

   // Generates the linearly spaced float vector ( 2.1, 3.2, 4.3, 5.4 )
   blaze::DynamicVector<float,columnVector> b;
   b = generate( 4UL, []( size_t index ){ return 2.1F + 1.1F*index; } );

   // Generates the logarithmically spaced double vector ( 1.0, 10.0, 100.0, 1000.0 )
   blaze::DynamicVector<double,columnVector> c;
   c = generate<columnVector>( 4UL, []( size_t index ){ return blaze::exp10( 1.0 + 1.0*index ); } );

   // Generates the vector of integer vectors ( ( 1, 2 ), ( 2, 3 ), ( 3, 4 ), ( 4, 5 ) )
   using VT = StaticVector<int,2UL>;
   blaze::DynamicVector<VT,rowVector> d;
   d = generate<rowVector>( 4UL, []( size_t index ) { return evaluate( VT{ 1, 2 } + index ); } );
   \endcode
*/
template< bool TF = defaultTransposeFlag  // Transpose flag
        , typename OP >                   // Type of the custom unary operation
inline decltype(auto) generate( size_t size, OP op )
{
   BLAZE_FUNCTION_TRACE;

   using ET = RemoveCVRef_t< decltype( std::declval<OP>()( std::declval<size_t>() ) ) >;
   using VT = DynamicVector<ET,TF>;

   return generate<VT>( size, std::move( op ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Generates a new dense vector filled with linearly spaced elements.
// \ingroup dense_vector
//
// \param size The size/dimension of the given vector.
// \param start The value of the first element.
// \param end The value of the last element.
// \return The newly generated dense vector.
//
// The \a linspace() function returns a dense vector filled with linearly spaced elements from
// \a start to \a end. By default, the returned vector is a column vector, but this setting can
// be changed via the \c BLAZE_DEFAULT_TRANSPOSE_FLAG switch. Alternatively it is possible to
// specify the transpose flag explicitly.\n
// The following example demonstrates the use of the \a linspace() function:

   \code
   using blaze::linspace;
   using blaze::columnVector;
   using blaze::rowVector;

   // Generates the linearly spaced integer vector ( 2, 3, 4, 5, 6 )
   blaze::DynamicVector<int,columnVector> a;
   a = linspace( 5UL, 2, 6 );

   // Generates the linearly spaced integer vector ( 6, 5, 4, 3, 2 )
   blaze::DynamicVector<int,columnVector> b;
   b = linspace<columnVector>( 5UL, 6, 2 );

   // Generates the linearly spaced float vector ( 2.1, 3.2, 4.3, 5.4 )
   blaze::DynamicVector<float,rowVector> c;
   c = linspace<rowVector>( 4UL, 2.1F, 5.4F );
   \endcode
*/
template< bool TF = defaultTransposeFlag  // Transpose flag
        , typename T >                    // Type of the delimiters
inline decltype(auto) linspace( size_t size, T start, T end )
{
   BLAZE_FUNCTION_TRACE;

   using BT = UnderlyingBuiltin_t<T>;
   using ET = If_t< IsFloatingPoint_v<BT>, BT, double >;

   ET divisor{ 1.0 };

   if( size > 2UL ) {
      divisor = static_cast<ET>( size - 1UL );
   }
   else {
      start = end;
   }

   auto delta( evaluate( ( end - start ) / divisor ) );

   return generate<TF>( size, [ start=std::move(start), delta=std::move(delta) ]( size_t index ) {
      return evaluate( start + index*delta );
   } );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Generates a new dense vector filled with logarithmically spaced elements.
// \ingroup dense_vector
//
// \param size The size/dimension of the given vector.
// \param start The value of the first element.
// \param end The value of the last element.
// \return The newly generated dense vector.
//
// The \a logspace() function returns a dense vector filled with logarithmically spaced elements
// from \a start to \a end. By default, the returned vector is a column vector, but this setting
// can be changed via the \c BLAZE_DEFAULT_TRANSPOSE_FLAG switch. Alternatively it is possible to
// specify the transpose flag explicitly.\n
// The following example demonstrates the use of the \a logspace() function:

   \code
   using blaze::logspace;
   using blaze::columnVector;
   using blaze::rowVector;

   // Generates the logarithmically spaced double vector ( 1, 10, 100, 1000 )
   blaze::DynamicVector<int,columnVector> a;
   a = logspace( 4UL, 0, 3 );

   // Generates the logarithmically spaced double vector ( 1000.0, 100.0, 10.0, 1.0 )
   blaze::DynamicVector<double,rowVector> b;
   b = logspace<rowVector>( 4UL, 3.0, 0.0 );
   \endcode
*/
template< bool TF = defaultTransposeFlag  // Transpose flag
        , typename T >                    // Type of the delimiters
inline decltype(auto) logspace( size_t size, T start, T end )
{
   BLAZE_FUNCTION_TRACE;

   using BT = UnderlyingBuiltin_t<T>;
   using ET = If_t< IsFloatingPoint_v<BT>, BT, double >;

   ET divisor{ 1.0 };

   if( size > 2UL ) {
      divisor = static_cast<ET>( size - 1UL );
   }
   else {
      start = end;
   }

   auto delta( evaluate( ( end - start ) / divisor ) );

   return generate<TF>( size, [ start=std::move(start), delta=std::move(delta) ]( size_t index ) {
      return evaluate( exp10( start + index*delta ) );
   } );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS (SUBVECTOR)
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given dense vector generator expression.
// \ingroup dense_vector
//
// \param expr The dense vector generator expression containing the subvector.
// \param args The optional subvector arguments.
// \return View on the specified subvector of the dense vector generator expression.
// \exception std::invalid_argument Invalid subvector specification.
//
// This function returns an expression representing the specified subvector of the given dense
// vector generator expression.
*/
template< AlignmentFlag AF    // Alignment flag
        , size_t I            // Index of the first subvector element
        , size_t N            // Size of the subvector
        , typename VT         // Type of the dense vector
        , typename OP         // Type of the custom unary operation
        , bool TF             // Transpose flag
        , typename... RSAs >  // Optional arguments
inline decltype(auto) subvector( const DVecGenExpr<VT,OP,TF>& expr, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( args... );

   constexpr bool isChecked( !Contains_v< TypeList<RSAs...>, Unchecked > );

   if( isChecked ) {
      if( I + N > expr.size() ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid subvector specification" );
      }
   }
   else {
      BLAZE_USER_ASSERT( I + N <= expr.size(), "Invalid subvector specification" );
   }

   return generate( N, [op=expr.operation()]( size_t i ) {
      return op( i+I );
   } );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given dense vector generator expression.
// \ingroup dense_vector
//
// \param expr The dense vector generator expression containing the subvector.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \param args The optional subvector arguments.
// \return View on the specified subvector of the dense vector generator expression.
// \exception std::invalid_argument Invalid subvector specification.
//
// This function returns an expression representing the specified subvector of the given dense
// vector generator expression.
*/
template< AlignmentFlag AF    // Alignment flag
        , typename VT         // Type of the dense vector
        , typename OP         // Type of the custom unary operation
        , bool TF             // Transpose flag
        , typename... RSAs >  // Optional arguments
inline decltype(auto) subvector( const DVecGenExpr<VT,OP,TF>& expr, size_t index, size_t size, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( args... );

   constexpr bool isChecked( !Contains_v< TypeList<RSAs...>, Unchecked > );

   if( isChecked ) {
      if( index + size > expr.size() ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid subvector specification" );
      }
   }
   else {
      BLAZE_USER_ASSERT( index + size <= expr.size(), "Invalid subvector specification" );
   }

   return generate( size, [op=expr.operation(),index]( size_t i ) {
      return op( i+index );
   } );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS (ELEMENTS)
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements on a dense vector generator expression.
// \ingroup dense_vector
//
// \param expr The given dense vector generator expression.
// \param args The optional element arguments.
// \return View on the specified selection of elements on the dense vector generator expression.
// \exception std::invalid_argument Invalid elements specification.
//
// This function returns an expression representing the specified selection of elements on the
// given dense vector generator expression.
*/
template< size_t I            // First element index
        , size_t... Is        // Remaining element indices
        , typename VT         // Type of the dense vector
        , typename OP         // Type of the custom unary operation
        , bool TF             // Transpose flag
        , typename... REAs >  // Optional element arguments
inline decltype(auto) elements( const DVecGenExpr<VT,OP,TF>& expr, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( args... );

   constexpr bool isChecked( !Contains_v< TypeList<REAs...>, Unchecked > );

   if( isChecked ) {
      static constexpr size_t indices[] = { I, Is... };
      for( size_t i=0UL; i<sizeof...(Is)+1UL; ++i ) {
         if( expr.size() <= indices[i] ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid elements specification" );
         }
      }
   }

   return generate( sizeof...(Is)+1UL, [op=expr.operation()]( size_t i ) {
      static constexpr size_t indices[] = { I, Is... };
      return op( indices[i] );
   } );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements on a dense vector generator expression.
// \ingroup dense_vector
//
// \param expr The given dense vector generator expression.
// \param indices The container of element indices.
// \param n The total number of indices.
// \param args The optional element arguments.
// \return View on the specified selection of elements on the dense vector generator expression.
// \exception std::invalid_argument Invalid elements specification.
//
// This function returns an expression representing the specified selection of elements on the
// given dense vector generator expression.
*/
template< typename VT         // Type of the dense vector
        , typename OP         // Type of the custom unary operation
        , bool TF             // Transpose flag
        , typename T          // Type of the element indices
        , typename... REAs >  // Optional element arguments
inline decltype(auto) elements( const DVecGenExpr<VT,OP,TF>& expr, T* indices, size_t n, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( args... );

   constexpr bool isChecked( !Contains_v< TypeList<REAs...>, Unchecked > );

   if( isChecked ) {
      for( size_t i=0UL; i<n; ++i ) {
         if( expr.size() <= size_t( indices[i] ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid elements specification" );
         }
      }
   }

   SmallArray<size_t,128UL> newIndices( indices, indices+n );

   return generate( n, [op=expr.operation(),newIndices]( size_t i ) {
      return op( newIndices[i] );
   } );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements on a dense vector generator expression.
// \ingroup dense_vector
//
// \param expr The given dense vector generator expression.
// \param p Callable producing the indices.
// \param n The total number of indices.
// \param args The optional element arguments.
// \return View on the specified selection of elements on the dense vector generator expression.
// \exception std::invalid_argument Invalid elements specification.
//
// This function returns an expression representing the specified selection of elements on the
// given dense vector generator expression.
*/
template< typename VT         // Type of the dense vector
        , typename OP         // Type of the custom unary operation
        , bool TF             // Transpose flag
        , typename P          // Type of the index producer
        , typename... REAs >  // Optional element arguments
inline decltype(auto) elements( const DVecGenExpr<VT,OP,TF>& expr, P p, size_t n, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( args... );

   constexpr bool isChecked( !Contains_v< TypeList<REAs...>, Unchecked > );

   if( isChecked ) {
      for( size_t i=0UL; i<n; ++i ) {
         if( expr.size() <= size_t( p(i) ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid elements specification" );
         }
      }
   }

   return generate( n, [op=expr.operation(),p]( size_t i ) {
      return op( p( i ) );
   } );
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
