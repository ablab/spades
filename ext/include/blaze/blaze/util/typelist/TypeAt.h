//=================================================================================================
/*!
//  \file blaze/util/typelist/TypeAt.h
//  \brief Header file for the TypeAt class template
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

#ifndef _BLAZE_UTIL_TYPELIST_TYPEAT_H_
#define _BLAZE_UTIL_TYPELIST_TYPEAT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/InvalidType.h>
#include <blaze/util/typelist/TypeList.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Indexing a type list.
// \ingroup typelist
//
// The TypeAt class can be used to access a type list at a specified position to query the
// according type. In order to index a type list, the TypeAt class has to be instantiated
// for a particular type list and an index value. The indexed type is available via the
// member type definition \a Type. The following example gives an impression of the use
// of the TypeAt class:

   \code
   using Floats = blaze::TypeList< float, double, long double >;  // Defining a new type list
   using Index0 = blaze::TypeAt< Floats, 0 >::Type;               // Indexing of the type list at index 0
   \endcode

// \note The access index is zero based!
*/
template< typename TL     // Type of the type list
        , size_t Index >  // Type list access index
struct TypeAt;
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the TypeAt class for an index of 0.
// \ingroup typelist
*/
template< typename T        // Type at the head of the type list
        , typename... Ts >  // Types of the tail of the type list
struct TypeAt< TypeList<T,Ts...>, 0UL >
{
   using Type = T;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Spezialization of the TypeAt class for empty type lists.
// \ingroup typelist
*/
template< size_t Index >  // Type list access index
struct TypeAt< TypeList<>, Index >
{
   using Type = INVALID_TYPE;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the TypeAt class for a general index.
// \ingroup typelist
*/
template< typename T      // Type of the head of the type list
        , typename... Ts  // Types of the tail of the type list
        , size_t Index >  // Type list access index
struct TypeAt< TypeList<T,Ts...>, Index >
{
   using Type = typename TypeAt< TypeList<Ts...>, Index-1UL >::Type;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the TypeAt class template.
// \ingroup type_traits
//
// The TypeAt_t alias declaration provides a convenient shortcut to access the nested \a Type
// of the TypeAt class template. For instance, given the type list \a TL and the index \a Index
// the following two type definitions are identical:

   \code
   using Type1 = typename TypeAt<TL,Index>::Type;
   using Type2 = TypeAt_t<TL,Index>;
   \endcode
*/
template< typename TL     // Type of the type list
        , size_t Index >  // Type list access index
using TypeAt_t = typename TypeAt<TL,Index>::Type;
//*************************************************************************************************

} // namespace blaze

#endif
