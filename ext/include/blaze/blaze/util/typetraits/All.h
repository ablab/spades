//=================================================================================================
/*!
//  \file blaze/util/typetraits/All.h
//  \brief Header file for the All type trait
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

#ifndef _BLAZE_UTIL_TYPETRAITS_ALL_H_
#define _BLAZE_UTIL_TYPETRAITS_ALL_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/IntegralConstant.h>
#include <blaze/util/mpl/And.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time type check.
// \ingroup type_traits
//
// This type trait determines whether the given type trait \a TypeTrait evaluates to \a true for
// all given types \a Ts. If the expression

   \code
   And_t< TypeTrait<Ts>... >::value
   \endcode

// evaluates to \a true, the \a value member constant is set to \a true, the nested type definition
// \a Type is \a TrueType, and the class derives from \a TrueType. Otherwise \a value is set to
// \a false, \a Type is \a FalseType, and the class derives from \a FalseType. Examples:

   \code
   blaze::All< IsIntegral, int, short, long >::value      // Evaluates to 'true'
   blaze::All< IsPointer, int*, float* >::Type            // Results in TrueType
   blaze::All< IsCharacter, char, signed char, wchar_t >  // Is derived from TrueType
   blaze::All< IsIntegral, int, float, double >::value    // Evaluates to 'false'
   blaze::All< IsPointer, int*, float& >::Type            // Results in FalseType
   blaze::All< IsCharacter, char, signed int, wchar_t >   // Is derived from FalseType
   \endcode
*/
template< template< typename > class TypeTrait  // Type trait to be evaluated on all operands
        , typename T1                           // Type of the first mandatory operand
        , typename T2                           // Type of the second mandatory operand
        , typename... Ts >                      // Types of the optional operands
struct All
   : public BoolConstant< And_t< TypeTrait<T1>, TypeTrait<T2>, TypeTrait<Ts>... >::value >
{};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary variable template for the All type trait.
// \ingroup type_traits
//
// The All_v variable template provides a convenient shortcut to access the nested \a value
// of the All class template. For instance, given the type trait \a TypeTrait and the two
// types \a T1 and \a T2 the following two statements are identical:

   \code
   constexpr bool value1 = blaze::All<TypeTrait,T1,T2>::value;
   constexpr bool value2 = blaze::All_v<TypeTrait,T1,T2>;
   \endcode
*/
template< template< typename > class TypeTrait  // Type trait to be evaluated on all operands
        , typename T1                           // Type of the first mandatory operand
        , typename T2                           // Type of the second mandatory operand
        , typename... Ts >                      // Types of the optional operands
constexpr bool All_v = All<TypeTrait,T1,T2,Ts...>::value;
//*************************************************************************************************

} // namespace blaze

#endif
