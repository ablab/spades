//=================================================================================================
/*!
//  \file blaze/math/typetraits/IsPaddingEnabled.h
//  \brief Header file for the IsPaddingEnabled type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_ISPADDINGENABLED_H_
#define _BLAZE_MATH_TYPETRAITS_ISPADDINGENABLED_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/IntegralConstant.h>
#include <blaze/util/typetraits/AlwaysFalse.h>
#include <blaze/util/typetraits/HasMember.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief First auxiliary helper struct for the IsPaddingEnabled type trait.
// \ingroup math_type_traits
*/
BLAZE_CREATE_HAS_DATA_OR_FUNCTION_MEMBER_TYPE_TRAIT( IsPaddingEnabledHelper1, paddingEnabled );
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Second auxiliary helper struct for the IsPaddingEnabled type trait.
// \ingroup math_type_traits
*/
template< typename T >
struct IsPaddingEnabledHelper2 {
   static constexpr bool test( bool (*fnc)() ) { return fnc(); }
   static constexpr bool test( bool b ) { return b; }
   static constexpr bool value = test( T::paddingEnabled );
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compile time check for data types.
// \ingroup math_type_traits
//
// This type trait queries the nested \a paddingEnabled member of the given data type \a T, which
// indicates whether the type provides support for padding (i.e. can properly deal with zeros).
// If the type is supporting padding, the \a value member constant is set to \a true, the nested
// type definition \a Type is \a TrueType, and the class derives from \a TrueType. Otherwise
// \a value is set to \a false, \a Type is \a FalseType, and the class derives from \a FalseType.
// Examples:

   \code
   struct A { static constexpr bool paddingEnabled = true; };
   struct B { static constexpr bool paddingEnabled() { return true; } };

   struct C {};
   struct D { static constexpr bool paddingEnabled = false; };
   struct E { static constexpr bool paddingEnabled() { return false; } };

   blaze::IsPaddingEnabled< A >::value  // Evaluates to 1
   blaze::IsPaddingEnabled< A >::Type   // Results in TrueType
   blaze::IsPaddingEnabled< B >         // Is derived from TrueType
   blaze::IsPaddingEnabled< C >::value  // Evaluates to 0
   blaze::IsPaddingEnabled< D >::Type   // Results in FalseType
   blaze::IsPaddingEnabled< E >         // Is derived from FalseType
   \endcode
*/
template< typename T >
struct IsPaddingEnabled
   : public BoolConstant< If_t< IsPaddingEnabledHelper1_v<T>
                              , IsPaddingEnabledHelper2<T>
                              , AlwaysFalse<T> >::value >
{};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary variable template for the IsPaddingEnabled type trait.
// \ingroup math_type_traits
//
// The IsPaddingEnabled_v variable template provides a convenient shortcut to access the nested
// \a value of the IsPaddingEnabled class template. For instance, given the type \a T the
// following two statements are identical:

   \code
   constexpr bool value1 = blaze::IsPaddingEnabled<T>::value;
   constexpr bool value2 = blaze::IsPaddingEnabled_v<T>;
   \endcode
*/
template< typename T >
constexpr bool IsPaddingEnabled_v = IsPaddingEnabled<T>::value;
//*************************************************************************************************

} // namespace blaze

#endif
