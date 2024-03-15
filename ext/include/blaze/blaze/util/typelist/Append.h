//=================================================================================================
/*!
//  \file blaze/util/typelist/Append.h
//  \brief Header file for the Append class template
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

#ifndef _BLAZE_UTIL_TYPELIST_APPEND_H_
#define _BLAZE_UTIL_TYPELIST_APPEND_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/typelist/TypeList.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Appending a type to a type list.
// \ingroup typelist
//
// The Append class can be used to append the data type \a T to a type list \a TL. In order to
// append a data type, the Append class has to be instantiated for a particular type list and
// another type. The following example gives an impression of the use of the Append class:

   \code
   using Tmp    = blaze::TypeList< float, double >;      // Defining a temporary type list
   using Floats = blaze::Append<Tmp,long double>::Type;  // Type list contains all floating point data types
   \endcode
*/
template< typename TL   // Type of the type list
        , typename T >  // The type to be appended to the type list
struct Append;
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Spezialization of the Append class for appending a single type.
// \ingroup typelist
*/
template< typename... Ts  // Types of the type list
        , typename T >    // The type to be appended to the type list
struct Append< TypeList<Ts...>, T >
{
   using Type = TypeList<Ts...,T>;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Spezialization of the Append class for appending a type list.
// \ingroup typelist
*/
template< typename... Ts1    // Type of the type list
        , typename... Ts2 >  // The types to be appended to the type list
struct Append< TypeList<Ts1...>, TypeList<Ts2...> >
{
   using Type = TypeList<Ts1...,Ts2...>;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the Append class template.
// \ingroup type_traits
//
// The Append_t alias declaration provides a convenient shortcut to access the nested \a Type
// of the Append class template. For instance, given the type list \a TL and the type \a T the
// following two type definitions are identical:

   \code
   using Type1 = typename Append<TL,T>::Type;
   using Type2 = Append_t<TL,T>;
   \endcode
*/
template< typename TL   // Type of the type list
        , typename T >  // The type to be appended to the type list
using Append_t = typename Append<TL,T>::Type;
//*************************************************************************************************

} // namespace blaze

#endif
