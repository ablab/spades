//=================================================================================================
/*!
//  \file blaze/util/typetraits/HasGreaterThan.h
//  \brief Header file for the HasGreaterThan type trait
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

#ifndef _BLAZE_UTIL_TYPETRAITS_HASGREATERTHAN_H_
#define _BLAZE_UTIL_TYPETRAITS_HASGREATERTHAN_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/util/typetraits/IsDetected.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary type alias for the HasGreaterThan type trait.
// \ingroup type_traits
*/
template< typename T1, typename T2 >
using GreaterThan_t = decltype( std::declval<T1>() > std::declval<T2>() );
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compile time check for the availability of a greater-than operation between two data types.
// \ingroup type_traits
//
// This type trait determines whether the given data types \a T1 and \a T2 can be used in a
// greater-than operation. If the operation is available, the \a value member constant is set
// to \a true, the nested type definition \a Type is \a TrueType, and the class derives from
// \a TrueType. Otherwise \a value is set to \a false, \a Type is \a FalseType, and the class
// derives from \a FalseType.

   \code
   blaze::HasGreaterThan< int, int >::value                       // Evaluates to 1
   blaze::HasGreaterThan< const std::string, std::string >::Type  // Results in TrueType
   blaze::HasGreaterThan< volatile int*, int* >                   // Is derived from TrueType
   blaze::HasGreaterThan< int, blaze::complex<float> >::value     // Evaluates to 0
   blaze::HasGreaterThan< std::string, int >::Type                // Results in FalseType
   blaze::HasGreaterThan< int*, std::string* >                    // Is derived from FalseType
   \endcode
*/
template< typename T1, typename T2 >
using HasGreaterThan = IsDetected< GreaterThan_t, T1, T2 >;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary variable template for the HasGreaterThan type trait.
// \ingroup type_traits
//
// The HasGreaterThan_v variable template provides a convenient shortcut to access the nested
// \a value of the HasGreaterThan class template. For instance, given the types \a T1 and \a T2
// the following two statements are identical:

   \code
   constexpr bool value1 = blaze::HasGreaterThan<T1,T2>::value;
   constexpr bool value2 = blaze::HasGreaterThan_v<T1,T2>;
   \endcode
*/
template< typename T1, typename T2 >
constexpr bool HasGreaterThan_v = HasGreaterThan<T1,T2>::value;
//*************************************************************************************************

} // namespace blaze

#endif
