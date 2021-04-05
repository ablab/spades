//=================================================================================================
/*!
//  \file blaze/math/typetraits/UnderlyingElement.h
//  \brief Header file for the UnderlyingElement type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_UNDERLYINGELEMENT_H_
#define _BLAZE_MATH_TYPETRAITS_UNDERLYINGELEMENT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/Complex.h>
#include <blaze/util/typetraits/RemoveCV.h>
#include <blaze/util/typetraits/Void.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename, typename = void > struct UnderlyingElementHelper1;
template< typename, typename = void > struct UnderlyingElementHelper2;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Evaluation of the element type of a given data type.
// \ingroup math_type_traits
//
// This type trait evaluates the underlying element type of the given data type \a T. If the given
// type provides a nested type \a ElementType, this type is reported as underlying element type
// type via the nested type \a Type. Else if the type provides a nested \a value_type, this type
// is reported as underlying element type. Else the given type itself reported as the underlying
// element type. Examples:

   \code
   using Type1 = double;                                    // Built-in data type
   using Type2 = complex<float>;                            // Complex data type
   using Type3 = std::vector<short>;                        // std::vector with built-in element type
   using Type4 = StaticVector<int,3UL>;                     // Vector with built-in element type
   using Type5 = CompressedMatrix< DynamicVector<float> >;  // Matrix with vector element type

   blaze::UnderlyingElement< Type1 >::Type  // corresponds to double
   blaze::UnderlyingElement< Type2 >::Type  // corresponds to float
   blaze::UnderlyingElement< Type3 >::Type  // corresponds to short
   blaze::UnderlyingElement< Type4 >::Type  // corresponds to int
   blaze::UnderlyingElement< Type5 >::Type  // corresponds to DynamicVector<float>
   \endcode

// Note that it is possible to add support for other data types that have an underlying element
// type but do neither provide a nested \a ElementType nor \a value_type type by specializing
// the UnderlyingElement class template.
*/
template< typename T >
struct UnderlyingElement
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   using Type = typename UnderlyingElementHelper1< RemoveCV_t<T> >::Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the UnderlyingElement type trait.
// \ingroup math_type_traits
//
// The UnderlyingElement_t alias declaration provides a convenient shortcut to access the
// nested \a Type of the UnderlyingElement class template. For instance, given the type \a T
// the following two type definitions are identical:

   \code
   using Type1 = typename blaze::UnderlyingElement<T>::Type;
   using Type2 = blaze::UnderlyingElement_t<T>;
   \endcode
*/
template< typename T >
using UnderlyingElement_t = typename UnderlyingElement<T>::Type;
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief First auxiliary helper struct for the UnderlyingElement type trait.
// \ingroup math_type_traits
*/
template< typename T, typename >
struct UnderlyingElementHelper1
{
   using Type = typename UnderlyingElementHelper2<T>::Type;
};

template< typename T >
struct UnderlyingElementHelper1< complex<T> >
{
   using Type = T;
};

template< typename T >
struct UnderlyingElementHelper1< T, Void_t< typename T::ElementType > >
{
   using Type = typename T::ElementType;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Second auxiliary helper struct for the UnderlyingElement type trait.
// \ingroup math_type_traits
*/
template< typename T, typename >
struct UnderlyingElementHelper2
{
   using Type = T;
};

template< typename T >
struct UnderlyingElementHelper2< T, Void_t< typename T::value_type > >
{
   using Type = typename T::value_type;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
