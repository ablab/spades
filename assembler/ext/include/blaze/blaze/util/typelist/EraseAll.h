//=================================================================================================
/*!
//  \file blaze/util/typelist/EraseAll.h
//  \brief Header file for the EraseAll class template
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

#ifndef _BLAZE_UTIL_TYPELIST_ERASEALL_H_
#define _BLAZE_UTIL_TYPELIST_ERASEALL_H_


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
/*!\brief Erasing all occurrences of a type from a type list.
// \ingroup typelist
//
// The EraseAll class can be used to erase all occurrences of data type \a Type from a type list
// \a TList. In order to erase all occurrences of a data type, the EraseAll class has to be
// instantiated for a particular type list and another type. The following example gives an
// impression of the use of the EraseAll class:

   \code
   // Defining a temporary type list containing the type int twice
   using Tmp = blaze::TypeList< float, int, double, int >;

   // Erasing the all occurrences of int from the type list
   using NoInt = blaze::EraseAll<Tmp,int>::Type;
   \endcode
*/
template< typename TL   // Type of the type list
        , typename T >  // The type to be erased from the type list
struct EraseAll;
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Spezialization of the EraseAll class for empty type lists.
// \ingroup typelist
*/
template< typename T >  // The type to be erased from the type list
struct EraseAll< TypeList<>, T >
{
   using Type = TypeList<>;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Spezialization of the EraseAll class for erasing an occurrence of T.
// \ingroup typelist
*/
template< typename T        // The type to be erased from the type list
        , typename... Ts >  // Type of the tail of the type list
struct EraseAll< TypeList<T,Ts...>, T >
{
   using Type = typename EraseAll< TypeList<Ts...>, T >::Type;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Spezialization of the EraseAll class for a general type list.
// \ingroup typelist
*/
template< typename U      // Type of the head of the type list
        , typename... Ts  // Types of the tail of the type list
        , typename T >    // The search type
struct EraseAll< TypeList<U,Ts...>, T >
{
   using Type = typename Append< TypeList<U>, typename EraseAll< TypeList<Ts...>, T >::Type >::Type;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the EraseAll class template.
// \ingroup type_traits
//
// The EraseAll_t alias declaration provides a convenient shortcut to access the nested \a Type
// of the EraseAll class template. For instance, given the type list \a TL and the type \a T the
// following two type definitions are identical:

   \code
   using Type1 = typename EraseAll<TL,T>::Type;
   using Type2 = EraseAll_t<TL,T>;
   \endcode
*/
template< typename TL   // Type of the type list
        , typename T >  // The type to be erased from the type list
using EraseAll_t = typename EraseAll<TL,T>::Type;
//*************************************************************************************************

} // namespace blaze

#endif
