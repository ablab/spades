//=================================================================================================
/*!
//  \file blaze/math/typetraits/UnderlyingBuiltin.h
//  \brief Header file for the UnderlyingBuiltin type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_UNDERLYINGBUILTIN_H_
#define _BLAZE_MATH_TYPETRAITS_UNDERLYINGBUILTIN_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/EnableIf.h>
#include <blaze/util/typetraits/IsSame.h>
#include <blaze/util/typetraits/RemoveCV.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename, typename = void > struct UnderlyingBuiltinHelper1;
template< typename, typename = void > struct UnderlyingBuiltinHelper2;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Evaluation of the underlying builtin element type of a given data type.
// \ingroup math_type_traits
//
// This type trait evaluates the underlying fundamental element type at the heart of the given
// data type \a T. For this purpose either a nested \a ElementType or a nested \a value_type will
// be used. Examples:

   \code
   using Type1 = double;                                    // Built-in data type
   using Type2 = complex<float>;                            // Complex data type
   using Type3 = std::vector<short>;                        // std::vector with built-in element type
   using Type4 = StaticVector<int,3UL>;                     // Vector with built-in element type
   using Type5 = CompressedVector< DynamicVector<float> >;  // Vector with vector element type

   blaze::UnderlyingBuiltin< Type1 >::Type  // corresponds to double
   blaze::UnderlyingBuiltin< Type2 >::Type  // corresponds to float
   blaze::UnderlyingBuiltin< Type3 >::Type  // corresponds to short
   blaze::UnderlyingBuiltin< Type4 >::Type  // corresponds to int
   blaze::UnderlyingBuiltin< Type5 >::Type  // corresponds to float
   \endcode

// Note that it is possible to add support for other data types that have an underlying
// fundamental element type but do neither provide a nested \a ElementType nor \a value_type type
// by specializing the UnderlyingBuiltin class template.
*/
template< typename T >
struct UnderlyingBuiltin
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   using Type = typename UnderlyingBuiltinHelper1< RemoveCV_t<T> >::Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the UnderlyingBuiltin type trait.
// \ingroup math_type_traits
//
// The UnderlyingBuiltin_t alias declaration provides a convenient shortcut to access the
// nested \a Type of the UnderlyingBuiltin class template. For instance, given the type \a T
// the following two type definitions are identical:

   \code
   using Type1 = typename blaze::UnderlyingBuiltin<T>::Type;
   using Type2 = blaze::UnderlyingBuiltin_t<T>;
   \endcode
*/
template< typename T >
using UnderlyingBuiltin_t = typename UnderlyingBuiltin<T>::Type;
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief First auxiliary helper struct for the UnderlyingBuiltin type trait.
// \ingroup math_type_traits
*/
template< typename T, typename >
struct UnderlyingBuiltinHelper1
{
   using Type = typename UnderlyingBuiltinHelper2<T>::Type;
};

template< typename T >
struct UnderlyingBuiltinHelper1< T, EnableIf_t< !IsSame_v< T, typename T::ElementType > > >
{
   using Type = typename UnderlyingBuiltinHelper1< typename T::ElementType >::Type;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Second auxiliary helper struct for the UnderlyingBuiltin type trait.
// \ingroup math_type_traits
*/
template< typename T, typename >
struct UnderlyingBuiltinHelper2
{
   using Type = T;
};

template< typename T >
struct UnderlyingBuiltinHelper2< T, EnableIf_t< !IsSame_v< T, typename T::value_type > > >
{
   using Type = typename UnderlyingBuiltinHelper1< typename T::value_type >::Type;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
