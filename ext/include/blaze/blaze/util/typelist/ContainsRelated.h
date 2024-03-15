//=================================================================================================
/*!
//  \file blaze/util/typelist/ContainsRelated.h
//  \brief Header file for the ContainsRelated class template
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

#ifndef _BLAZE_UTIL_TYPELIST_CONTAINSRELATED_H_
#define _BLAZE_UTIL_TYPELIST_CONTAINSRELATED_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/IntegralConstant.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/typelist/TypeList.h>
#include <blaze/util/typetraits/IsConvertible.h>


namespace blaze {

//=================================================================================================
//
//  CLASS TEMPLATE
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Searching a type list.
// \ingroup typelist
//
// The ContainsRelated class can be used to search the type list for a type related to \a Type.
// In contrast to the Contains class, the ContainsRelated class only searches for a type the
// given data type \a Type can be converted to. In case a related type is found in the type list,
// the \a value member enumeration is set to \a true, else it is set to \a false. In order to
// check whether a related type is contained in the type list, the ContainsRelated class has to
// be instantiated for a particular type list and another type. The following example gives an
// impression of the use of the ContainsRelated class:

   \code
   class A {};
   class B : public A {};
   class C {};
   class D {};

   // Defining a new type list
   using Types = blaze::TypeList< A, C >;

   // Searching for the type A in the type list
   constexpr bool a = blaze::ContainsRelated< Types, A >::value;  // Evaluates to 1, type A is found

   // Searching for the derived type B in the type list
   constexpr bool b = blaze::ContainsRelated< Types, B >::value;  // Evaluates to 1, base type A is found

   // Searching for the type C in the type list
   constexpr bool c = blaze::ContainsRelated< Types, D >::value;  // Evaluates to 0, no related type found
   \endcode
*/
template< typename TL   // Type of the type list
        , typename T >  // The search type
struct ContainsRelated;
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Spezialization of the ContainsRelated class for empty type lists.
// \ingroup typelist
*/
template< typename T >  // The search type
struct ContainsRelated< TypeList<>, T >
   : public FalseType
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Spezialization of the ContainsRelated class for a general type list.
// \ingroup typelist
*/
template< typename U      // Type of the head of the type list
        , typename... Ts  // Types of the tail of the type list
        , typename T >    // The search type
struct ContainsRelated< TypeList<U,Ts...>, T >
   : public If_t< IsConvertible_v<T,U>
                , TrueType
                , ContainsRelated< TypeList<Ts...>, T > >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary variable template for the ContainsRelated type trait.
// \ingroup type_traits
//
// The ContainsRelated_v variable template provides a convenient shortcut to access the nested
// \a value of the ContainsRelated class template. For instance, given the type list \a TL and
// the type \a T the following two statements are identical:

   \code
   constexpr bool value1 = ContainsRelated<TL,T>::value;
   constexpr bool value2 = ContainsRelated_v<TL,T>;
   \endcode
*/
template< typename TL   // Type of the type list
        , typename T >  // The search type
constexpr bool ContainsRelated_v = ContainsRelated<TL,T>::value;
//*************************************************************************************************

} // namespace blaze

#endif
