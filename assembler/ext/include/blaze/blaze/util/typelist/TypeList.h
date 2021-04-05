//=================================================================================================
/*!
//  \file blaze/util/typelist/TypeList.h
//  \brief Header file for the TypeList class template
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

#ifndef _BLAZE_UTIL_TYPELIST_TYPELIST_H_
#define _BLAZE_UTIL_TYPELIST_TYPELIST_H_


namespace blaze {

//=================================================================================================
//
//  CLASS TYPELIST
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup typelist Type lists
// \ingroup util
//
// Type lists provide the functionality to create lists of data types. In constrast to lists of
// data values (as for instance the std::list class template), type lists are created at compile
// time, not at run time. The following example demonstrates, how type lists are created and
// manipulated:

   \code
   // Creating a type list consisting of two fundamental floating point data types
   using Tmp = blaze::TypeList< float, double >;

   // Appending a type to the type list
   using Floats = blaze::Append< Tmp, long double >::Type;  // Type list contains all floating point data types

   // Calculating the length of the type list (at compile time!)
   const int length = blaze::Length< Floats >::value;  // Value evaluates to 3

   // Accessing a specific type of the type list via indexing
   using Index0 = blaze::TypeAt< Floats, 0 >::Type;

   // Searching the type list for a specific type
   constexpr bool index1 = blaze::Contains< Floats, double >::value;  // Value evaluates to 1
   constexpr bool index2 = blaze::Contains< Floats, int    >::value;  // Value evaluates to 0

   // Estimating the index of a specific type in the type list
   constexpr bool index3 = blaze::IndexOf< Floats, double >::value;   // Value evaluates to 1
   constexpr bool index4 = blaze::IndexOf< Floats, int    >::value;   // Value evaluates to 3

   // Erasing the first occurrence of float from the type list
   using NoFloat = blaze::Erase< Floats, float >::Type;

   // Removing all duplicates from the type list
   using NoDuplicates = blaze::Unique< Floats >::Type;
   \endcode
*/
/*!\brief Implementation of a type list.
// \ingroup typelist
//
// The TypeList class template represents a list of data types of arbitrary size. The following
// example gives an impression how type lists are used and manipulated:

   \code
   // Creating a type list consisting of two fundamental floating point data types
   using Tmp = blaze::TypeList< float, double >;

   // Appending a type to the type list
   using Floats = blaze::Append< Tmp, long double >::Type;  // Type list contains all floating point data types

   // Calculating the length of the type list (at compile time!)
   const int length = blaze::Length< Floats >::value;  // Value evaluates to 3

   // Accessing a specific type of the type list via indexing
   using Index0 = blaze::TypeAt< Floats, 0 >::Type;

   // Searching the type list for a specific type
   constexpr bool index1 = blaze::Contains< Floats, double >::value;  // Value evaluates to 1
   constexpr bool index2 = blaze::Contains< Floats, int    >::value;  // Value evaluates to 0

   // Estimating the index of a specific type in the type list
   constexpr bool index3 = blaze::IndexOf< Floats, double >::value;   // Value evaluates to 1
   constexpr bool index4 = blaze::IndexOf< Floats, int    >::value;   // Value evaluates to 3

   // Erasing the first occurrence of float from the type list
   using NoFloat = blaze::Erase< Floats, float >::Type;

   // Removing all duplicates from the type list
   using NoDuplicates = blaze::Unique< Floats >::Type;
   \endcode
*/
template< typename... Ts >
struct TypeList
{};
//*************************************************************************************************

} // namespace blaze

#endif
