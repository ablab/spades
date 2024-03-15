//=================================================================================================
/*!
//  \file blaze/util/typelist/Erase.h
//  \brief Header file for the Erase class template
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

#ifndef _BLAZE_UTIL_TYPELIST_ERASE_H_
#define _BLAZE_UTIL_TYPELIST_ERASE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/typelist/Append.h>
#include <blaze/util/typelist/TypeList.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Erasing the first occurrence of a type from a type list.
// \ingroup typelist
//
// The Erase class can be used to erase the first occurrence of data type \a T from a type
// list \a TL. In order to erase the first occurrence of a data type, the Erase class has to
// be instantiated for a particular type list and another type. The following example gives
// an impression of the use of the Erase class:

   \code
   // Defining a temporary type list containing the type int twice
   using Tmp = blaze::TypeList< float, int, double, int >;

   // Erasing the first occurrence of int from the type list
   using SingleInt = blaze::Erase<Tmp,int>::Type;
   \endcode
*/
template< typename TL   // Type of the type list
        , typename T >  // The type to be erased from the type list
struct Erase;
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Spezialization of the Erase class for empty type lists.
// \ingroup typelist
*/
template< typename T >  // The type to be erased from the type list
struct Erase< TypeList<>, T >
{
   using Type = TypeList<>;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Spezialization of the Erase class for erasing the first occurrence of T.
// \ingroup typelist
*/
template< typename T        // The type to be erased from the type list
        , typename... Ts >  // Types of the tail of the type list
struct Erase< TypeList<T,Ts...>, T >
{
   using Type = TypeList<Ts...>;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Spezialization of the Erase class for a general type list.
// \ingroup typelist
*/
template< typename U      // Type of the head of the type list
        , typename... Ts  // Types of the tail of the type list
        , typename T >    // The search type
struct Erase< TypeList<U,Ts...>, T >
{
   using Type = typename Append< TypeList<U>, typename Erase< TypeList<Ts...>, T >::Type >::Type;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the Erase class template.
// \ingroup type_traits
//
// The Erase_t alias declaration provides a convenient shortcut to access the nested \a Type
// of the Erase class template. For instance, given the type list \a TL and the type \a T the
// following two type definitions are identical:

   \code
   using Type1 = typename Erase<TL,T>::Type;
   using Type2 = Erase_t<TL,T>;
   \endcode
*/
template< typename TL   // Type of the type list
        , typename T >  // The type to be appended to the type list
using Erase_t = typename Erase<TL,T>::Type;
//*************************************************************************************************

} // namespace blaze

#endif
