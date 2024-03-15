//=================================================================================================
/*!
//  \file blaze/math/typetraits/MaxSize.h
//  \brief Header file for the MaxSize type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_MAXSIZE_H_
#define _BLAZE_MATH_TYPETRAITS_MAXSIZE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/IntegralConstant.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/Void.h>


namespace blaze {

//=================================================================================================
//
//  ::blaze NAMESPACE FORWARD DECLARATIONS
//
//=================================================================================================

template< typename, size_t, typename = void > struct MaxSizeHelper1;
template< typename, size_t, typename > struct MaxSizeHelper2;




//=================================================================================================
//
//  TYPE DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Default size of the MaxSize type trait.
// \ingroup math_type_traits
*/
constexpr ptrdiff_t DefaultMaxSize_v = -1L;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Type representation of the default size of the MaxSize type trait.
// \ingroup math_type_traits
*/
using DefaultMaxSize = Ptrdiff_t<DefaultMaxSize_v>;
//*************************************************************************************************




//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time evaluation of the maximum size of vectors and matrices.
// \ingroup math_type_traits
//
// The MaxSize type trait evaluates the maximum size of a particular dimension of the given
// vector or matrix type at compile time. In case the given type \a T is a vector or matrix type
// with a fixed maximum size (e.g. StaticVector, HybridVector, StaticMatrix, or HybridMatrix)
// and \a N is a valid dimension, the \a value member constant is set to the according size.
// In all other cases, \a value is set to -1.

   \code
   using blaze::StaticVector;
   using blaze::HybridMatrix;
   using blaze::DynamicVector;

   blaze::MaxSize< StaticVector<int,3UL>, 0UL >::value      // Evaluates to 3
   blaze::MaxSize< HybridMatrix<int,2UL,4UL>, 0UL >::value  // Evaluates to 2 (the number of rows)
   blaze::MaxSize< HybridMatrix<int,2UL,4UL>, 1UL >::value  // Evaluates to 4 (the number of columns)
   blaze::MaxSize< StaticVector<int,3UL>, 1UL >::value      // Evaluates to -1; 1 is not a valid vector dimension!
   blaze::MaxSize< DynamicVector<int>, 0UL >::value         // Evaluates to -1; Maximum size not fixed at compile time!
   blaze::MaxSize< int, 0UL >::value                        // Evaluates to -1
   \endcode
*/
template< typename T, size_t N >
struct MaxSize
   : public MaxSizeHelper1<T,N>
{};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the MaxSize type trait for const types.
// \ingroup math_type_traits
*/
template< typename T, size_t N >
struct MaxSize< const T, N >
   : public MaxSize<T,N>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the MaxSize type trait for volatile types.
// \ingroup math_type_traits
*/
template< typename T, size_t N >
struct MaxSize< volatile T, N >
   : public MaxSize<T,N>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the MaxSize type trait for cv qualified types.
// \ingroup math_type_traits
*/
template< typename T, size_t N >
struct MaxSize< const volatile T, N >
   : public MaxSize<T,N>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary variable template for the MaxSize type trait.
// \ingroup math_type_traits
//
// The MaxSize_v variable template provides a convenient shortcut to access the nested \a value
// of the MaxSize class template. For instance, given the type \a T and the dimension \a N the
// following two statements are identical:

   \code
   constexpr size_t value1 = blaze::MaxSize<T,N>::value;
   constexpr size_t value2 = blaze::MaxSize_v<T,N>;
   \endcode
*/
template< typename T, size_t N >
constexpr ptrdiff_t MaxSize_v = MaxSize<T,N>::value;
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief First auxiliary helper struct for the MaxSize type trait.
// \ingroup math_type_traits
*/
template< typename T, size_t N, typename >
struct MaxSizeHelper1
   : public DefaultMaxSize
{};

template< typename T, size_t N >
struct MaxSizeHelper1< T, N, Void_t< typename T::ResultType > >
   : public MaxSizeHelper2<T,N,typename T::ResultType>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Second auxiliary helper struct for the MaxSize type trait.
// \ingroup math_type_traits
*/
template< typename T, size_t N, typename U >
struct MaxSizeHelper2
   : public MaxSize<U,N>
{};

template< typename T, size_t N >
struct MaxSizeHelper2<T,N,T>
   : public DefaultMaxSize
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
