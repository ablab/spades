//=================================================================================================
/*!
//  \file blaze/util/typelist/IndexOf.h
//  \brief Header file for the IndexOf class template
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

#ifndef _BLAZE_UTIL_TYPELIST_INDEXOF_H_
#define _BLAZE_UTIL_TYPELIST_INDEXOF_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/IntegralConstant.h>
#include <blaze/util/typelist/TypeList.h>


namespace blaze {

//=================================================================================================
//
//  TYPE LIST SEARCH
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Searching a type list.
// \ingroup typelist
//
// The IndexOf class can be used to search the type list for a particular type \a Type. In
// contrast to the Contains and the ContainsRelated classes, the IndexOf class evaluates the
// index of the given type in the type list. In case the type is contained in the type list,
// the \a value member represents the index of the queried type. Otherwise the \a value member
// is set to the length of the type list. In order to search for a type, the IndexOf class has
// to be instantiated for a particular type list and a search type. The following example gives
// an impression of the use of the IndexOf class:

   \code
   using Floats = blaze::TypeList< float, double, long double >;     // Defining a new type list
   constexpr bool index1 = blaze::IndexOf< Floats, double >::value;  // Value evaluates to 1
   constexpr bool index2 = blaze::IndexOf< Floats, int    >::value;  // Value evaluates to -1
   \endcode
*/
template< typename TL   // Type of the type list
        , typename T >  // The search type
struct IndexOf;
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Spezialization of the IndexOf class for empty type lists.
// \ingroup typelist
*/
template< typename T >  // The search type
struct IndexOf< TypeList<>, T >
   : public Size_t<1UL>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Spezialization of the IndexOf class for a successful search.
// \ingroup typelist
*/
template< typename T        // The search type
        , typename... Ts >  // Types of the tail of the type list
struct IndexOf< TypeList<T,Ts...>, T >
   : public Size_t<0UL>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Spezialization of the IndexOf class for a general type list.
// \ingroup typelist
*/
template< typename U      // Type of the head of the type list
        , typename... Ts  // Types of the tail of the type list
        , typename T >    // The search type
struct IndexOf< TypeList<U,Ts...>, T >
   : public Size_t< 1UL + IndexOf< TypeList<Ts...>, T >::value >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary variable template for the IndexOf type trait.
// \ingroup type_traits
//
// The IndexOf_v variable template provides a convenient shortcut to access the nested \a value
// of the IndexOf class template. For instance, given the type list \a TL and the type \a T the
// following two statements are identical:

   \code
   constexpr size_t value1 = IndexOf<TL,T>::value;
   constexpr size_t value2 = IndexOf_v<TL,T>;
   \endcode
*/
template< typename TL   // Type of the type list
        , typename T >  // The search type
constexpr size_t IndexOf_v = IndexOf<TL,T>::value;
//*************************************************************************************************

} // namespace blaze

#endif
