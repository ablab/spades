//=================================================================================================
/*!
//  \file blaze/math/expressions/Vector.h
//  \brief Header file for the Vector CRTP base class
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

#ifndef _BLAZE_MATH_EXPRESSIONS_VECTOR_H_
#define _BLAZE_MATH_EXPRESSIONS_VECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Exception.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsShrinkable.h>
#include <blaze/system/Inline.h>
#include <blaze/system/MacroDisable.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsSame.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup vector Vectors
// \ingroup math
*/
/*!\brief Base class for N-dimensional vectors.
// \ingroup vector
//
// The Vector class is a base class for all arbitrarily sized (N-dimensional) dense and sparse
// vector classes within the Blaze library. It provides an abstraction from the actual type of
// the vector, but enables a conversion back to this type via the 'Curiously Recurring Template
// Pattern' (CRTP).
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag
class Vector
{
 public:
   //**Type definitions****************************************************************************
   using VectorType = VT;  //!< Type of the vector.
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   static constexpr bool transposeFlag = TF;  //!< Transpose flag of the vector.
   //**********************************************************************************************

   //**Conversion operators************************************************************************
   /*!\name Conversion operators */
   //@{
   [[deprecated]] BLAZE_ALWAYS_INLINE constexpr VT&       operator~()       noexcept;
   [[deprecated]] BLAZE_ALWAYS_INLINE constexpr const VT& operator~() const noexcept;

   constexpr VT&       operator*()       noexcept;
   constexpr const VT& operator*() const noexcept;
   //@}
   //**********************************************************************************************

 protected:
   //**Special member functions********************************************************************
   /*!\name Special member functions */
   //@{
   Vector() = default;
   Vector( const Vector& ) = default;
   Vector( Vector&& ) = default;
   ~Vector() = default;
   Vector& operator=( const Vector& ) = default;
   Vector& operator=( Vector&& ) = default;
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CONVERSION OPERATIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief CRTP-based conversion operation for non-constant vectors.
//
// \return Mutable reference of the actual type of the vector.
//
// This operator performs the CRTP-based type-safe downcast to the actual type \a VT of the
// vector. It will return a mutable reference to the actual type \a VT.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag
[[deprecated]] BLAZE_ALWAYS_INLINE constexpr VT& Vector<VT,TF>::operator~() noexcept
{
   return static_cast<VT&>( *this );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief CRTP-based conversion operation for constant vectors.
//
// \return Constant reference of the actual type of the vector.
//
// This operator performs the CRTP-based type-safe downcast to the actual type \a VT of the
// vector. It will return a constant reference to the actual type \a VT.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag
[[deprecated]] BLAZE_ALWAYS_INLINE constexpr const VT& Vector<VT,TF>::operator~() const noexcept
{
   return static_cast<const VT&>( *this );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief CRTP-based conversion operation for non-constant vectors.
//
// \return Mutable reference of the actual type of the vector.
//
// This operator performs the CRTP-based type-safe downcast to the actual type \a VT of the
// vector. It will return a mutable reference to the actual type \a VT.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag
BLAZE_ALWAYS_INLINE constexpr VT& Vector<VT,TF>::operator*() noexcept
{
   return static_cast<VT&>( *this );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief CRTP-based conversion operation for constant vectors.
//
// \return Const reference of the actual type of the vector.
//
// This operator performs the CRTP-based type-safe downcast to the actual type \a VT of the
// vector. It will return a constant reference to the actual type \a VT.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag
BLAZE_ALWAYS_INLINE constexpr const VT& Vector<VT,TF>::operator*() const noexcept
{
   return static_cast<const VT&>( *this );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Vector global functions */
//@{
template< typename VT, bool TF >
VT& crtp_cast( Vector<VT,TF>& vector );

template< typename VT, bool TF >
const VT& crtp_cast( const Vector<VT,TF>& vector );

template< typename VT, bool TF >
typename VT::Iterator begin( Vector<VT,TF>& vector );

template< typename VT, bool TF >
typename VT::ConstIterator begin( const Vector<VT,TF>& vector );

template< typename VT, bool TF >
typename VT::ConstIterator cbegin( const Vector<VT,TF>& vector );

template< typename VT, bool TF >
typename VT::Iterator end( Vector<VT,TF>& vector );

template< typename VT, bool TF >
typename VT::ConstIterator end( const Vector<VT,TF>& vector );

template< typename VT, bool TF >
typename VT::ConstIterator cend( const Vector<VT,TF>& vector );

template< typename VT, bool TF >
constexpr size_t size( const Vector<VT,TF>& vector ) noexcept;

template< typename VT, bool TF >
size_t capacity( const Vector<VT,TF>& vector ) noexcept;

template< typename VT, bool TF >
size_t nonZeros( const Vector<VT,TF>& vector );

template< typename VT, bool TF >
void resize( Vector<VT,TF>& vector, size_t n, bool preserve=true );

template< typename VT, bool TF >
void shrinkToFit( Vector<VT,TF>& vector );

template< typename VT, bool TF >
typename VT::ResultType evaluate( const Vector<VT,TF>& vector );

template< bool B, typename VT, bool TF >
decltype(auto) evaluateIf( const Vector<VT,TF>& vector );

template< typename VT, bool TF >
constexpr bool isEmpty( const Vector<VT,TF>& vector ) noexcept;

template< typename VT1, bool TF1, typename VT2, bool TF2 >
bool isSame( const Vector<VT1,TF1>& a, const Vector<VT2,TF2>& b ) noexcept;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief CRTP-based conversion operation for non-constant vectors.
//
// \param vector The vector to be downcast.
// \return Mutable reference of the actual type of the vector.
//
// This operator performs the CRTP-based type-safe downcast to the actual type \a VT of the
// vector. It will return a mutable reference to the actual type \a VT.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
BLAZE_ALWAYS_INLINE VT& crtp_cast( Vector<VT,TF>& vector )
{
   return *vector;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief CRTP-based conversion operation for constant vectors.
//
// \param vector The vector to be downcast.
// \return Const reference of the actual type of the vector.
//
// This operator performs the CRTP-based type-safe downcast to the actual type \a VT of the
// vector. It will return a constant reference to the actual type \a VT.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
BLAZE_ALWAYS_INLINE const VT& crtp_cast( const Vector<VT,TF>& vector )
{
   return *vector;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the given vector.
// \ingroup vector
//
// \param vector The given dense or sparse vector.
// \return Iterator to the first element of the given vector.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
BLAZE_ALWAYS_INLINE typename VT::Iterator begin( Vector<VT,TF>& vector )
{
   return (*vector).begin();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the given vector.
// \ingroup vector
//
// \param vector The given dense or sparse vector.
// \return Iterator to the first element of the given vector.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
BLAZE_ALWAYS_INLINE typename VT::ConstIterator begin( const Vector<VT,TF>& vector )
{
   return (*vector).begin();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the given vector.
// \ingroup vector
//
// \param vector The given dense or sparse vector.
// \return Iterator to the first element of the given vector.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
BLAZE_ALWAYS_INLINE typename VT::ConstIterator cbegin( const Vector<VT,TF>& vector )
{
   return (*vector).begin();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the given vector.
// \ingroup vector
//
// \param vector The given dense or sparse vector.
// \return Iterator just past the last element of the given vector.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
BLAZE_ALWAYS_INLINE typename VT::Iterator end( Vector<VT,TF>& vector )
{
   return (*vector).end();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the given vector.
// \ingroup vector
//
// \param vector The given dense or sparse vector.
// \return Iterator just past the last element of the given vector.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
BLAZE_ALWAYS_INLINE typename VT::ConstIterator end( const Vector<VT,TF>& vector )
{
   return (*vector).end();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the given vector.
// \ingroup vector
//
// \param vector The given dense or sparse vector.
// \return Iterator just past the last element of the given vector.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
BLAZE_ALWAYS_INLINE typename VT::ConstIterator cend( const Vector<VT,TF>& vector )
{
   return (*vector).end();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current size/dimension of the vector.
// \ingroup vector
//
// \param vector The given vector.
// \return The size of the vector.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
BLAZE_ALWAYS_INLINE constexpr size_t size( const Vector<VT,TF>& vector ) noexcept
{
   return (*vector).size();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the maximum capacity of the vector.
// \ingroup vector
//
// \param vector The given vector.
// \return The capacity of the vector.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
BLAZE_ALWAYS_INLINE size_t capacity( const Vector<VT,TF>& vector ) noexcept
{
   return (*vector).capacity();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of non-zero elements in the vector.
// \ingroup vector
//
// \param vector The given vector.
// \return The number of non-zero elements in the vector.
//
// Note that the number of non-zero elements is always less than or equal to the current size
// of the vector.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
BLAZE_ALWAYS_INLINE size_t nonZeros( const Vector<VT,TF>& vector )
{
   return (*vector).nonZeros();
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c resize() function for non-resizable vectors.
// \ingroup vector
//
// \param vector The given vector to be resized.
// \param n The new size of the vector.
// \return void
// \exception std::invalid_argument Vector cannot be resized.
//
// This function tries to change the number of rows and columns of a non-resizable vector. Since
// the vector cannot be resized, in case the specified size is not identical to the current size
// of the vector, a \a std::invalid_argument exception is thrown.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
BLAZE_ALWAYS_INLINE auto resize_backend( Vector<VT,TF>& vector, size_t n, bool preserve )
   -> DisableIf_t< IsResizable_v<VT> >
{
   MAYBE_UNUSED( preserve );

   if( (*vector).size() != n ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector cannot be resized" );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c resize() function for resizable vectors.
// \ingroup vector
//
// \param vector The given vector to be resized.
// \param n The new size of the vector.
// \param preserve \a true if the old values of the vector should be preserved, \a false if not.
// \return void
//
// This function changes the size of the given resizable vector.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
BLAZE_ALWAYS_INLINE auto resize_backend( Vector<VT,TF>& vector, size_t n, bool preserve )
   -> EnableIf_t< IsResizable_v<VT> >
{
   (*vector).resize( n, preserve );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Changing the size of the vector.
// \ingroup vector
//
// \param vector The given vector to be resized.
// \param n The new size of the vector.
// \param preserve \a true if the old values of the vector should be preserved, \a false if not.
// \return void
// \exception std::invalid_argument Vector cannot be resized.
//
// This function provides a unified interface to resize dense and sparse vectors. In contrast
// to the \c resize() member function, which is only available on resizable vector types, this
// function can be used on both resizable and non-resizable vectors. In case the given vector
// type \a VT is resizable (i.e. provides a \c resize() function), the type-specific \c resize()
// member function is called. Depending on the type \a VT, this may result in the allocation of
// new dynamic memory and the invalidation of existing views (subvectors, ...). In case \a VT is
// non-resizable (i.e. does not provide a \c resize() function) and if the specified size is not
// identical to the current size of the vector, a \a std::invalid_argument exception is thrown.

   \code
   blaze::DynamicVector<int> a( 3UL );
   resize( a, 5UL );  // OK: regular resize operation

   blaze::StaticVector<int,3UL> b;
   resize( b, 3UL );  // OK: No resize necessary
   resize( b, 5UL );  // Error: Vector cannot be resized!
   \endcode
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
BLAZE_ALWAYS_INLINE void resize( Vector<VT,TF>& vector, size_t n, bool preserve )
{
   resize_backend( vector, n, preserve );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c shrinkToFit() function for non-shrinkable vectors.
// \ingroup vector
//
// \param vector The given vector to be shrunk.
// \return void
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
BLAZE_ALWAYS_INLINE auto shrinkToFit_backend( Vector<VT,TF>& vector )
   -> DisableIf_t< IsShrinkable_v<VT> >
{
   MAYBE_UNUSED( vector );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c shrinkToFit() function for shrinkable vectors.
// \ingroup vector
//
// \param vector The given vector to be shrunk.
// \return void
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
BLAZE_ALWAYS_INLINE auto shrinkToFit_backend( Vector<VT,TF>& vector )
   -> EnableIf_t< IsShrinkable_v<VT> >
{
   (*vector).shrinkToFit();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Requesting the removal of unused capacity.
// \ingroup vector
//
// \param vector The given vector to be shrunk.
// \return void
//
// This function tries to minimize the capacity of the vector by removing unused capacity.
// Please note that in case of a shrinkable vector, due to padding the capacity might not be
// reduced exactly to the size of the vector. Please also note that in case a reallocation
// occurs, all iterators (including end() iterators), all pointers and references to elements
// of this vector are invalidated. In case of an unshrinkable vector the function has no effect.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
BLAZE_ALWAYS_INLINE void shrinkToFit( Vector<VT,TF>& vector )
{
   shrinkToFit_backend( vector );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Evaluates the given vector expression.
// \ingroup vector
//
// \param vector The vector to be evaluated.
// \return The result of the evaluated vector expression.
//
// This function forces an evaluation of the given vector expression and enables an automatic
// deduction of the correct result type of an operation. The following code example demonstrates
// its intended use for the multiplication of a dense and a sparse vector:

   \code
   using blaze::DynamicVector;
   using blaze::CompressedVector;

   blaze::DynamicVector<double> a;
   blaze::CompressedVector<double> b;
   // ... Resizing and initialization

   auto c = evaluate( a * b );
   \endcode

// In this scenario, the \a evaluate() function assists in deducing the exact result type of
// the operation via the 'auto' keyword. Please note that if \a evaluate() is used in this
// way, no temporary vector is created and no copy operation is performed. Instead, the result
// is directly written to the target vector due to the return value optimization (RVO). However,
// if \a evaluate() is used in combination with an explicit target type, a temporary will be
// created and a copy operation will be performed if the used type differs from the type
// returned from the function:

   \code
   CompressedVector<double> d( a * b );  // No temporary & no copy operation
   DynamicVector<double> e( a * b );     // Temporary & copy operation
   d = evaluate( a * b );                // Temporary & copy operation
   \endcode

// Sometimes it might be desirable to explicitly evaluate a sub-expression within a larger
// expression. However, please note that \a evaluate() is not intended to be used for this
// purpose. This task is more elegantly and efficiently handled by the \a eval() function:

   \code
   blaze::DynamicVector<double> a, b, c, d;

   d = a + evaluate( b * c );  // Unnecessary creation of a temporary vector
   d = a + eval( b * c );      // No creation of a temporary vector
   \endcode

// In contrast to the \a evaluate() function, \a eval() can take the complete expression into
// account and therefore can guarantee the most efficient way to evaluate it.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
inline typename VT::ResultType evaluate( const Vector<VT,TF>& vector )
{
   typename VT::ResultType tmp( *vector );
   return tmp;
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conditional evaluation the given vector expression.
// \ingroup vector
//
// \param vector The vector to be evaluated.
// \return The result of the evaluated vector expression.
//
// This function does not evaluate the given vector expression and returns a reference to the
// vector expression.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
inline decltype(auto) evaluateIf( FalseType, const Vector<VT,TF>& vector )
{
   return *vector;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conditional evaluation the given vector expression.
// \ingroup vector
//
// \param vector The vector to be evaluated.
// \return The result of the evaluated vector expression.
//
// This function evaluates the given vector expression by means of the evaluate() function.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
inline decltype(auto) evaluateIf( TrueType, const Vector<VT,TF>& vector )
{
   return evaluate( *vector );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conditional evaluation of the given vector expression.
// \ingroup vector
//
// \param vector The vector to be evaluated.
// \return The result of the evaluated vector expression.
//
// In case the given compile time condition evaluates to \a true, this function evaluates the
// the given vector expression by means of the evaluate() function. Otherwise the function returns
// a reference to the given vector expression.
*/
template< bool B       // Compile time condition
        , typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
inline decltype(auto) evaluateIf( const Vector<VT,TF>& vector )
{
   return evaluateIf( BoolConstant<B>{}, *vector );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given vector is empty.
// \ingroup vector
//
// \param vector The vector to be checked.
// \return \a true if the vector is empty, \a false if not.
//
// This function checks if the total number of elements of the given vector is zero. If the
// total number of elements is zero the function returns \a true, otherwise it returns \a false.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
BLAZE_ALWAYS_INLINE constexpr bool isEmpty( const Vector<VT,TF>& vector ) noexcept
{
   return size( *vector ) == 0UL;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the two given vectors represent the same observable state.
// \ingroup vector
//
// \param a The first vector to be tested for its state.
// \param b The second vector to be tested for its state.
// \return \a true in case the two vectors share a state, \a false otherwise.
//
// The isSame function provides an abstract interface for testing if the two given vectors
// represent the same observable state. This happens for instance in case \c a and \c b refer
// to the same vector or in case \c a and \c b are aliases for the same vector. In case both
// vectors represent the same observable state, the function returns \a true, other it returns
// \a false.

   \code
   blaze::DynamicVector<int> vec1( 4UL );  // Setup of a 4-dimensional dynamic vector
   blaze::DynamicVector<int> vec2( 4UL );  // Setup of a second 4-dimensional dynamic vector

   auto sub1 = subvector( vec1, 0UL, 4UL );  // Subvector of vec1 for the entire range
   auto sub2 = subvector( vec1, 1UL, 2UL );  // Subvector of vec1 for the range [1..3]
   auto sub3 = subvector( vec1, 1UL, 2UL );  // Second subvector of vec1 for the range [1..3]

   isSame( vec1, vec1 );  // returns true since both objects refer to the same vector
   isSame( vec1, vec2 );  // returns false since vec1 and vec2 are two different vectors
   isSame( vec1, sub1 );  // returns true since sub1 represents the same observable state as vec1
   isSame( vec1, sub3 );  // returns false since sub3 only covers part of the range of vec1
   isSame( sub2, sub3 );  // returns true since sub1 and sub2 refer to exactly the same range of vec1
   isSame( sub1, sub3 );  // returns false since sub1 and sub3 refer to different ranges of vec1
   \endcode
*/
template< typename VT1  // Type of the left-hand side vector
        , bool TF1      // Transpose flag of the left-hand side vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
BLAZE_ALWAYS_INLINE bool isSame( const Vector<VT1,TF1>& a, const Vector<VT2,TF2>& b ) noexcept
{
   return ( IsSame_v<VT1,VT2> &&
            reinterpret_cast<const void*>( &a ) == reinterpret_cast<const void*>( &b ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a vector to a vector.
// \ingroup vector
//
// \param lhs The target left-hand side vector.
// \param rhs The right-hand side vector to be assigned.
// \return void
//
// This function implements the default assignment of a vector to another vector.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side vector
        , bool TF1      // Transpose flag of the left-hand side vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
BLAZE_ALWAYS_INLINE void assign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).size() == (*rhs).size(), "Invalid vector sizes" );
   (*lhs).assign( *rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a vector to a vector.
// \ingroup vector
//
// \param lhs The target left-hand side vector.
// \param rhs The right-hand side vector to be added.
// \return void
//
// This function implements the default addition assignment of a vector to a vector.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side vector
        , bool TF1      // Transpose flag of the left-hand side vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
BLAZE_ALWAYS_INLINE void addAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).size() == (*rhs).size(), "Invalid vector sizes" );
   (*lhs).addAssign( *rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a vector to a vector.
// \ingroup vector
//
// \param lhs The target left-hand side vector.
// \param rhs The right-hand side vector to be subtracted.
// \return void
//
// This function implements the default subtraction assignment of a vector to a vector.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side vector
        , bool TF1      // Transpose flag of the left-hand side vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
BLAZE_ALWAYS_INLINE void subAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).size() == (*rhs).size(), "Invalid vector sizes" );
   (*lhs).subAssign( *rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the multiplication assignment of a vector to a vector.
// \ingroup vector
//
// \param lhs The target left-hand side vector.
// \param rhs The right-hand side vector to be multiplied.
// \return void
//
// This function implements the default multiplication assignment of a vector to a vector.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side vector
        , bool TF1      // Transpose flag of the left-hand side vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
BLAZE_ALWAYS_INLINE void multAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).size() == (*rhs).size(), "Invalid vector sizes" );
   (*lhs).multAssign( *rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the division assignment of a vector to a vector.
// \ingroup vector
//
// \param lhs The target left-hand side vector.
// \param rhs The right-hand side vector divisor.
// \return void
//
// This function implements the default division assignment of a vector to a vector.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side vector
        , bool TF1      // Transpose flag of the left-hand side vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
BLAZE_ALWAYS_INLINE void divAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).size() == (*rhs).size(), "Invalid vector sizes" );
   (*lhs).divAssign( *rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by setting a single element of a vector.
// \ingroup vector
//
// \param vec The target vector.
// \param index The index of the element to be set.
// \param value The value to be set to the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT    // Type of the vector
        , bool TF        // Transpose flag
        , typename ET >  // Type of the element
BLAZE_ALWAYS_INLINE bool trySet( const Vector<VT,TF>& vec, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < (*vec).size(), "Invalid vector access index" );

   MAYBE_UNUSED( vec, index, value );

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by setting a range of elements of a vector.
// \ingroup vector
//
// \param vec The target vector.
// \param index The index of the first element of the range to be set.
// \param size The number of elements of the range to be set.
// \param value The value to be set to the range of elements.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT    // Type of the vector
        , bool TF        // Transpose flag
        , typename ET >  // Type of the element
BLAZE_ALWAYS_INLINE bool
   trySet( const Vector<VT,TF>& vec, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (*vec).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*vec).size(), "Invalid range size" );

   MAYBE_UNUSED( vec, index, size, value );

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by adding to a single element of a vector.
// \ingroup vector
//
// \param vec The target vector.
// \param index The index of the element to be modified.
// \param value The value to be added to the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT    // Type of the vector
        , bool TF        // Transpose flag
        , typename ET >  // Type of the element
BLAZE_ALWAYS_INLINE bool tryAdd( const Vector<VT,TF>& vec, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < (*vec).size(), "Invalid vector access index" );

   MAYBE_UNUSED( vec, index, value );

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by adding to a range of elements of a vector.
// \ingroup vector
//
// \param vec The target vector.
// \param index The index of the first element of the range to be modified.
// \param size The number of elements of the range to be modified.
// \param value The value to be added to the range of elements.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT    // Type of the vector
        , bool TF        // Transpose flag
        , typename ET >  // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryAdd( const Vector<VT,TF>& vec, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (*vec).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*vec).size(), "Invalid range size" );

   MAYBE_UNUSED( vec, index, size, value );

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by subtracting from a single element of a vector.
// \ingroup vector
//
// \param vec The target vector.
// \param index The index of the element to be modified.
// \param value The value to be subtracted from the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT    // Type of the vector
        , bool TF        // Transpose flag
        , typename ET >  // Type of the element
BLAZE_ALWAYS_INLINE bool trySub( const Vector<VT,TF>& vec, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < (*vec).size(), "Invalid vector access index" );

   MAYBE_UNUSED( vec, index, value );

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by subtracting from a range of elements of a vector.
// \ingroup vector
//
// \param vec The target vector.
// \param index The index of the first element of the range to be modified.
// \param size The number of elements of the range to be modified.
// \param value The value to be subtracted from the range of elements.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT    // Type of the vector
        , bool TF        // Transpose flag
        , typename ET >  // Type of the element
BLAZE_ALWAYS_INLINE bool
   trySub( const Vector<VT,TF>& vec, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (*vec).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*vec).size(), "Invalid range size" );

   MAYBE_UNUSED( vec, index, size, value );

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a single element of a vector.
// \ingroup vector
//
// \param vec The target vector.
// \param index The index of the element to be modified.
// \param value The factor for the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT    // Type of the vector
        , bool TF        // Transpose flag
        , typename ET >  // Type of the element
BLAZE_ALWAYS_INLINE bool tryMult( const Vector<VT,TF>& vec, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < (*vec).size(), "Invalid vector access index" );

   MAYBE_UNUSED( vec, index, value );

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a range of elements of a vector.
// \ingroup vector
//
// \param vec The target vector.
// \param index The index of the first element of the range to be modified.
// \param size The number of elements of the range to be modified.
// \param value The factor for the elements.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT    // Type of the vector
        , bool TF        // Transpose flag
        , typename ET >  // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryMult( const Vector<VT,TF>& vec, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (*vec).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*vec).size(), "Invalid range size" );

   MAYBE_UNUSED( vec, index, size, value );

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a single element of a vector.
// \ingroup vector
//
// \param vec The target vector.
// \param index The index of the element to be modified.
// \param value The divisor for the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT    // Type of the vector
        , bool TF        // Transpose flag
        , typename ET >  // Type of the element
BLAZE_ALWAYS_INLINE bool tryDiv( const Vector<VT,TF>& vec, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < (*vec).size(), "Invalid vector access index" );

   MAYBE_UNUSED( vec, index, value );

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a range of elements of a vector.
// \ingroup vector
//
// \param vec The target vector.
// \param index The index of the first element of the range to be modified.
// \param size The number of elements of the range to be modified.
// \param value The divisor for the elements.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT    // Type of the vector
        , bool TF        // Transpose flag
        , typename ET >  // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryDiv( const Vector<VT,TF>& vec, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (*vec).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*vec).size(), "Invalid range size" );

   MAYBE_UNUSED( vec, index, size, value );

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by shifting a single element of a vector.
// \ingroup vector
//
// \param vec The target vector.
// \param index The index of the element to be modified.
// \param count The number of bits to shift the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag
BLAZE_ALWAYS_INLINE bool tryShift( const Vector<VT,TF>& vec, size_t index, int count )
{
   BLAZE_INTERNAL_ASSERT( index < (*vec).size(), "Invalid vector access index" );

   MAYBE_UNUSED( vec, index, count );

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by shifting a range of elements of a vector.
// \ingroup vector
//
// \param vec The target vector.
// \param index The index of the first element of the range to be modified.
// \param size The number of elements of the range to be modified.
// \param count The number of bits to shift the range of elements.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag
BLAZE_ALWAYS_INLINE bool
   tryShift( const Vector<VT,TF>& vec, size_t index, size_t size, int count )
{
   BLAZE_INTERNAL_ASSERT( index <= (*vec).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*vec).size(), "Invalid range size" );

   MAYBE_UNUSED( vec, index, size, count );

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise AND on a single element of a vector.
// \ingroup vector
//
// \param vec The target vector.
// \param index The index of the element to be modified.
// \param value The bit pattern to be used on the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT    // Type of the vector
        , bool TF        // Transpose flag
        , typename ET >  // Type of the element
BLAZE_ALWAYS_INLINE bool tryBitand( const Vector<VT,TF>& vec, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < (*vec).size(), "Invalid vector access index" );

   MAYBE_UNUSED( vec, index, value );

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise AND on a range of elements of a vector.
// \ingroup vector
//
// \param vec The target vector.
// \param index The index of the first element of the range to be modified.
// \param size The number of elements of the range to be modified.
// \param value The bit pattern to be used on the range of elements.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT    // Type of the vector
        , bool TF        // Transpose flag
        , typename ET >  // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryBitand( const Vector<VT,TF>& vec, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (*vec).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*vec).size(), "Invalid range size" );

   MAYBE_UNUSED( vec, index, size, value );

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise OR on a single element of a vector.
// \ingroup vector
//
// \param vec The target vector.
// \param index The index of the element to be modified.
// \param value The bit pattern to be used on the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT    // Type of the vector
        , bool TF        // Transpose flag
        , typename ET >  // Type of the element
BLAZE_ALWAYS_INLINE bool tryBitor( const Vector<VT,TF>& vec, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < (*vec).size(), "Invalid vector access index" );

   MAYBE_UNUSED( vec, index, value );

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise OR on a range of elements of a vector.
// \ingroup vector
//
// \param vec The target vector.
// \param index The index of the first element of the range to be modified.
// \param size The number of elements of the range to be modified.
// \param value The bit pattern to be used on the range of elements.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT    // Type of the vector
        , bool TF        // Transpose flag
        , typename ET >  // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryBitor( const Vector<VT,TF>& vec, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (*vec).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*vec).size(), "Invalid range size" );

   MAYBE_UNUSED( vec, index, size, value );

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise XOR on a single element of a vector.
// \ingroup vector
//
// \param vec The target vector.
// \param index The index of the element to be modified.
// \param value The bit pattern to be used on the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT    // Type of the vector
        , bool TF        // Transpose flag
        , typename ET >  // Type of the element
BLAZE_ALWAYS_INLINE bool tryBitxor( const Vector<VT,TF>& vec, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < (*vec).size(), "Invalid vector access index" );

   MAYBE_UNUSED( vec, index, value );

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise XOR on a range of elements of a vector.
// \ingroup vector
//
// \param vec The target vector.
// \param index The index of the first element of the range to be modified.
// \param size The number of elements of the range to be modified.
// \param value The bit pattern to be used on the range of elements.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT    // Type of the vector
        , bool TF        // Transpose flag
        , typename ET >  // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryBitxor( const Vector<VT,TF>& vec, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (*vec).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*vec).size(), "Invalid range size" );

   MAYBE_UNUSED( vec, index, size, value );

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a vector to a vector.
// \ingroup vector
//
// \param lhs The target left-hand side vector.
// \param rhs The right-hand side vector to be assigned.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side vector
        , bool TF1      // Transpose flag of the left-hand side vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
BLAZE_ALWAYS_INLINE bool tryAssign( const Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= (*lhs).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= (*lhs).size(), "Invalid vector size" );

   MAYBE_UNUSED( lhs, rhs, index );

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a vector to a vector.
// \ingroup vector
//
// \param lhs The target left-hand side vector.
// \param rhs The right-hand side vector to be added.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side vector
        , bool TF1      // Transpose flag of the left-hand side vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
BLAZE_ALWAYS_INLINE bool tryAddAssign( const Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= (*lhs).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= (*lhs).size(), "Invalid vector size" );

   MAYBE_UNUSED( lhs, rhs, index );

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a vector to a vector.
// \ingroup vector
//
// \param lhs The target left-hand side vector.
// \param rhs The right-hand side vector to be subtracted.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side vector
        , bool TF1      // Transpose flag of the left-hand side vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
BLAZE_ALWAYS_INLINE bool trySubAssign( const Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= (*lhs).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= (*lhs).size(), "Invalid vector size" );

   MAYBE_UNUSED( lhs, rhs, index );

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the multiplication assignment of a vector to a vector.
// \ingroup vector
//
// \param lhs The target left-hand side vector.
// \param rhs The right-hand side vector to be multiplied.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side vector
        , bool TF1      // Transpose flag of the left-hand side vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
BLAZE_ALWAYS_INLINE bool tryMultAssign( const Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= (*lhs).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= (*lhs).size(), "Invalid vector size" );

   MAYBE_UNUSED( lhs, rhs, index );

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the division assignment of a vector to a vector.
// \ingroup vector
//
// \param lhs The target left-hand side vector.
// \param rhs The right-hand side vector divisor.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side vector
        , bool TF1      // Transpose flag of the left-hand side vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
BLAZE_ALWAYS_INLINE bool tryDivAssign( const Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= (*lhs).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= (*lhs).size(), "Invalid vector size" );

   MAYBE_UNUSED( lhs, rhs, index );

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the shift assignment of a vector to a vector.
// \ingroup vector
//
// \param lhs The target left-hand side vector.
// \param rhs The right-hand side vector of bits to shift.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side vector
        , bool TF1      // Transpose flag of the left-hand side vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
BLAZE_ALWAYS_INLINE bool tryShiftAssign( const Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= (*lhs).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= (*lhs).size(), "Invalid vector size" );

   MAYBE_UNUSED( lhs, rhs, index );

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the bitwise AND assignment of a vector to a vector.
// \ingroup vector
//
// \param lhs The target left-hand side vector.
// \param rhs The right-hand side vector for the bitwise AND operation.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side vector
        , bool TF1      // Transpose flag of the left-hand side vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
BLAZE_ALWAYS_INLINE bool tryBitandAssign( const Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= (*lhs).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= (*lhs).size(), "Invalid vector size" );

   MAYBE_UNUSED( lhs, rhs, index );

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the bitwise OR assignment of a vector to a vector.
// \ingroup vector
//
// \param lhs The target left-hand side vector.
// \param rhs The right-hand side vector for the bitwise OR operation.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side vector
        , bool TF1      // Transpose flag of the left-hand side vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
BLAZE_ALWAYS_INLINE bool tryBitorAssign( const Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= (*lhs).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= (*lhs).size(), "Invalid vector size" );

   MAYBE_UNUSED( lhs, rhs, index );

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the bitwise XOR assignment of a vector to a vector.
// \ingroup vector
//
// \param lhs The target left-hand side vector.
// \param rhs The right-hand side vector for the bitwise XOR operation.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side vector
        , bool TF1      // Transpose flag of the left-hand side vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
BLAZE_ALWAYS_INLINE bool tryBitxorAssign( const Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= (*lhs).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= (*lhs).size(), "Invalid vector size" );

   MAYBE_UNUSED( lhs, rhs, index );

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given vector.
// \ingroup vector
//
// \param vector The vector to be derestricted.
// \return Reference to the vector without access restrictions.
//
// This function removes all restrictions on the data access to the given vector. It returns a
// reference to the vector that does provide the same interface but does not have any restrictions
// on the data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
BLAZE_ALWAYS_INLINE VT& derestrict( Vector<VT,TF>& vector )
{
   return *vector;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of the top-level view on the given vector.
// \ingroup vector
//
// \param vector The vector to be unviewed.
// \return Reference to the vector without view.
//
// This function removes the top-level view on the given vector and returns a reference to the
// unviewed vector.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
BLAZE_ALWAYS_INLINE VT& unview( Vector<VT,TF>& vector )
{
   return *vector;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of the top-level views on the given constant vector.
// \ingroup vector
//
// \param vector The constant vector to be unviewed.
// \return Reference to the vector without view.
//
// This function removes the top-level view on the given constant vector and returns a reference
// to the unviewed vector.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag of the vector
BLAZE_ALWAYS_INLINE const VT& unview( const Vector<VT,TF>& vector )
{
   return *vector;
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
