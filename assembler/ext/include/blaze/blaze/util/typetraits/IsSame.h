//=================================================================================================
/*!
//  \file blaze/util/typetraits/IsSame.h
//  \brief Header file for the IsSame and IsStrictlySame type traits
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

#ifndef _BLAZE_UTIL_TYPETRAITS_ISSAME_H_
#define _BLAZE_UTIL_TYPETRAITS_ISSAME_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/IntegralConstant.h>
#include <blaze/util/typetraits/RemoveCV.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time type relationship analysis.
// \ingroup type_traits
//
// This class tests if the two data types \a A and \a B are equal. For this type comparison,
// the cv-qualifiers of both data types are not ignored. If \a A and \a B are the same data
// type, then the \a value member constant is set to \a true, the nested type definition
// \a Type is \a TrueType, and the class derives from \a TrueType. Otherwise \a value is set
// to \a false, \a Type is \a FalseType, and the class derives from \a FalseType.

   \code
   blaze::IsStrictlySame<int,int>::value                   // Evaluates to 'true'
   blaze::IsStrictlySame<const double,const double>::Type  // Results in TrueType
   blaze::IsStrictlySame<volatile float,volatile float>    // Is derived from TrueType
   blaze::IsStrictlySame<char,wchar_t>::value              // Evaluates to 'false'
   blaze::IsStrictlySame<int,const int>::Type              // Results in FalseType
   blaze::IsStrictlySame<float,volatile float>             // Is derived from FalseType
   \endcode
*/
template< typename A, typename B >
struct IsStrictlySame
   : public FalseType
{};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Specialization of the IsStrictlySame class template for a single, matching data type.
template< typename T >
struct IsStrictlySame<T,T>
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary variable template for the IsStrictlySame type trait.
// \ingroup type_traits
//
// The IsStrictlySame_v variable template provides a convenient shortcut to access the nested
// \a value of the IsStrictlySame class template. For instance, given the types \a T1 and \a T2
// the following two statements are identical:

   \code
   constexpr bool value1 = blaze::IsStrictlySame<T1,T2>::value;
   constexpr bool value2 = blaze::IsStrictlySame_v<T1,T2>;
   \endcode
*/
template< typename A, typename B >
constexpr bool IsStrictlySame_v = IsStrictlySame<A,B>::value;
//*************************************************************************************************




//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Type relationship analysis.
// \ingroup type_traits
//
// This class tests if the two data types \a A and \a B are equal. For this type comparison,
// the cv-qualifiers of both data types are ignored. If \a A and \a B are the same data type
// (ignoring the cv-qualifiers), then the \a value member constant is set to \a true, the
// nested type definition \a Type is \a TrueType, and the class derives from \a TrueType.
// Otherwise \a value is set to \a false, \a Type is \a FalseType, and the class derives from
// \a FalseType.

   \code
   blaze::IsSame<int,int>::value               // Evaluates to 'true'
   blaze::IsSame<int,const int>::Type          // Results in TrueType
   blaze::IsSame<float,volatile float>         // Is derived from TrueType
   blaze::IsSame<char,wchar_t>::value          // Evaluates to 'false'
   blaze::IsSame<char,volatile float>::Type    // Results in FalseType
   blaze::IsSame<int,double>                   // Is derived from FalseType
   \endcode
*/
template< typename A, typename B >
struct IsSame
   : public BoolConstant< IsStrictlySame_v< RemoveCV_t<A>, RemoveCV_t<B> > >
{};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary variable template for the IsSame type trait.
// \ingroup type_traits
//
// The IsSame_v variable template provides a convenient shortcut to access the nested \a value
// of the IsSame class template. For instance, given the types \a T1 and \a T2 the following
// two statements are identical:

   \code
   constexpr bool value1 = blaze::IsSame<T1,T2>::value;
   constexpr bool value2 = blaze::IsSame_v<T1,T2>;
   \endcode
*/
template< typename A, typename B >
constexpr bool IsSame_v = IsSame<A,B>::value;
//*************************************************************************************************

} // namespace blaze

#endif
