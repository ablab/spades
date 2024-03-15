//=================================================================================================
/*!
//  \file blaze/math/typetraits/IsSIMDCombinable.h
//  \brief Header file for the IsSIMDCombinable type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_ISSIMDCOMBINABLE_H_
#define _BLAZE_MATH_TYPETRAITS_ISSIMDCOMBINABLE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/IntegralConstant.h>
#include <blaze/util/typetraits/IsIntegral.h>
#include <blaze/util/typetraits/IsNumeric.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the IsSIMDCombinable type trait.
// \ingroup math_type_traits
*/
template< typename T1
        , typename T2
        , typename... Ts >
struct IsSIMDCombinableHelper
   : public BoolConstant< IsSIMDCombinableHelper<T1,T2>::value &&
                          IsSIMDCombinableHelper<T2,Ts...>::value >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Specialization of the IsSIMDCombinableHelper class template for two matching types.
template< typename T >
struct IsSIMDCombinableHelper<T,T>
   : public BoolConstant< IsNumeric_v<T> >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Specialization of the IsSIMDCombinableHelper class template for two different types.
template< typename T1, typename T2 >
struct IsSIMDCombinableHelper<T1,T2>
   : public BoolConstant< IsNumeric_v<T1> && IsIntegral_v<T1> &&
                          IsNumeric_v<T2> && IsIntegral_v<T2> &&
                          sizeof(T1) == sizeof(T2) >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compile time check for data types.
// \ingroup math_type_traits
//
// This type trait tests whether or not the given types can be combined in the context of SIMD
// operations. By default, all numeric data types are considered combinable to themselves. Also,
// different integral types of the same size can be combined. If the types are determined to be
// SIMD combinable, the \a value member constant is set to \a true, the nested type definition
// \a Type is \a TrueType, and the class derives from \a TrueType. Otherwise \a value is set to
// \a false, \a Type is \a FalseType, and the class derives from \a FalseType.

   \code
   blaze::IsSIMDCombinable< int, unsigned int >::value        // Evaluates to 1
   blaze::IsSIMDCombinable< double, double >::Type            // Results in TrueType
   blaze::IsSIMDCombinable< complex<float>, complex<float> >  // Is derived from TrueType
   blaze::IsSIMDCombinable< int, float >::value               // Evaluates to 0
   blaze::IsSIMDCombinable< double, float >::Type             // Results in FalseType
   blaze::IsSIMDCombinable< complex<int>, complex<float> >    // Is derived from FalseType
   \endcode
*/
template< typename T1, typename T2, typename... Ts >
struct IsSIMDCombinable
   : public BoolConstant< IsSIMDCombinableHelper< T1, T2, Ts... >::value >
{};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary variable template for the IsSIMDCombinable type trait.
// \ingroup math_type_traits
//
// The IsSIMDCombinable_v variable template provides a convenient shortcut to access the nested
// \a value of the IsSIMDCombinable class template. For instance, given the types \a T1 and
// \a T2 the following two statements are identical:

   \code
   constexpr bool value1 = blaze::IsSIMDCombinable<T1,T2>::value;
   constexpr bool value2 = blaze::IsSIMDCombinable_v<T1,T2>;
   \endcode
*/
template< typename T1, typename T2, typename... Ts >
constexpr bool IsSIMDCombinable_v = IsSIMDCombinable<T1,T2,Ts...>::value;
//*************************************************************************************************

} // namespace blaze

#endif
