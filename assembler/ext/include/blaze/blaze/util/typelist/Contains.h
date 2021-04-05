//=================================================================================================
/*!
//  \file blaze/util/typelist/Contains.h
//  \brief Header file for the Contains class template
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

#ifndef _BLAZE_UTIL_TYPELIST_CONTAINS_H_
#define _BLAZE_UTIL_TYPELIST_CONTAINS_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/IntegralConstant.h>
#include <blaze/util/typelist/TypeList.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Searching a type list.
// \ingroup typelist
//
// The Contains class can be used to search the type list for a particular type \a Type. In
// contrast to the IndexOf class, the Contains class does not evaluate the index of the type but
// only checks whether or not the type is contained in the type list. Additionally, in contrast
// to the ContainsRelated class, the Contains class strictly searches for the given type \a Type
// and not for a related data type. In case the type is contained in the type list, the \a value
// member enumeration is set to \a true, else it is set to \a false. In order to check whether a
// type is part of a type list, the Contains class has to be instantiated for a particular type
// list and another type. The following example gives an impression of the use of the Contains
// class:

   \code
   using Floats = blaze::TypeList< float, double, long double >;      // Defining a new type list
   constexpr bool index1 = blaze::Contains< Floats, double >::value;  // Value evaluates to true
   constexpr bool index2 = blaze::Contains< Floats, int    >::value;  // Value evaluates to false
   \endcode
*/
template< typename TL   // Type of the type list
        , typename T >  // The search type
struct Contains;
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Spezialization of the Contains class for empty type lists.
// \ingroup typelist
*/
template< typename T >  // The search type
struct Contains< TypeList<>, T >
   : public FalseType
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Spezialization of the Contains class for a successful search.
// \ingroup typelist
*/
template< typename T        // The search type
        , typename... Ts >  // Types of the tail of the type list
struct Contains< TypeList<T,Ts...>, T >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Spezialization of the Contains class for a general type list.
// \ingroup typelist
*/
template< typename U      // Type of the head of the type list
        , typename... Ts  // Types of the tail of the type list
        , typename T >    // The search type
struct Contains< TypeList<U,Ts...>, T >
   : public Contains< TypeList<Ts...>, T >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary variable template for the Contains type trait.
// \ingroup type_traits
//
// The Contains_v variable template provides a convenient shortcut to access the nested \a value
// of the Contains class template. For instance, given the type list \a TL and the type \a T the
// following two statements are identical:

   \code
   constexpr bool value1 = Contains<TL,T>::value;
   constexpr bool value2 = Contains_v<TL,T>;
   \endcode
*/
template< typename TL   // Type of the type list
        , typename T >  // The search type
constexpr bool Contains_v = Contains<TL,T>::value;
//*************************************************************************************************

} // namespace blaze

#endif
