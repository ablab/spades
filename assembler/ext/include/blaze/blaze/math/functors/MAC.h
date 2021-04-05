//=================================================================================================
/*!
//  \file blaze/math/functors/MAC.h
//  \brief Header file for the MAC functor
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

#ifndef _BLAZE_MATH_FUNCTORS_MAC_H_
#define _BLAZE_MATH_FUNCTORS_MAC_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/math/constraints/SIMDPack.h>
#include <blaze/math/typetraits/HasSIMDAdd.h>
#include <blaze/math/typetraits/HasSIMDMult.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsUniform.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/typetraits/IsZero.h>
#include <blaze/math/typetraits/YieldsLower.h>
#include <blaze/math/typetraits/YieldsStrictlyLower.h>
#include <blaze/math/typetraits/YieldsStrictlyUpper.h>
#include <blaze/math/typetraits/YieldsUniform.h>
#include <blaze/math/typetraits/YieldsUpper.h>
#include <blaze/math/typetraits/YieldsZero.h>
#include <blaze/system/HostDevice.h>
#include <blaze/system/Inline.h>
#include <blaze/util/IntegralConstant.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Generic wrapper for the multiply-accumulate (\f$ (a*b)+c \f$; MAC) operation.
// \ingroup functors
*/
struct MAC
{
   //**********************************************************************************************
   /*!\brief Returns the result of the MAC operation for the given objects/values.
   //
   // \param a The first MAC operand.
   // \param b The second MAC operand.
   // \param c The third MAC operand
   // \return The result of the MAC operation for the given objects/values.
   */
   template< typename T1, typename T2, typename T3 >
   BLAZE_ALWAYS_INLINE BLAZE_DEVICE_CALLABLE decltype(auto)
      operator()( T1&& a, T2&& b, T3&& c ) const
   {
      return ( std::forward<T1>( a ) * std::forward<T2>( b ) ) + std::forward<T3>( c );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether SIMD is enabled for the specified data types \a T1, \a T2, and \a T3.
   //
   // \return \a true in case SIMD is enabled for the data types \a T1, \a T2, and \a T3, \a false if not.
   */
   template< typename T1, typename T2, typename T3 >
   static constexpr bool simdEnabled() { return HasSIMDMult_v<T1,T2> && HasSIMDAdd_v<T1,T3>; }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the operation supports padding, i.e. whether it can deal with zeros.
   //
   // \return \a true in case padding is supported, \a false if not.
   */
   static constexpr bool paddingEnabled() { return true; }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns the result of the MAC operation for the given SIMD vectors.
   //
   // \param a The first SIMD vector.
   // \param b The second SIMD vector.
   // \param c The third SIMD vector.
   // \return The result of the MAC operation for the given SIMD vectors.
   */
   template< typename T1, typename T2, typename T3 >
   BLAZE_ALWAYS_INLINE decltype(auto) load( const T1& a, const T2& b, const T3& c ) const
   {
      BLAZE_CONSTRAINT_MUST_BE_SIMD_PACK( T1 );
      BLAZE_CONSTRAINT_MUST_BE_SIMD_PACK( T2 );
      BLAZE_CONSTRAINT_MUST_BE_SIMD_PACK( T3 );
      return ( a * b ) + c;
   }
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  YIELDSUNIFORM SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2, typename T3 >
struct YieldsUniform<MAC,T1,T2,T3>
   : public BoolConstant< IsUniform_v<T1> && IsUniform_v<T2> && IsUniform_v<T3> >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  YIELDSLOWER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, typename MT3 >
struct YieldsLower<MAC,MT1,MT2,MT3>
   : public BoolConstant< IsLower_v<MT1> && IsLower_v<MT2> && IsLower_v<MT3> >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  YIELDSSTRICTLYLOWER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, typename MT3 >
struct YieldsStrictlyLower<MAC,MT1,MT2,MT3>
   : public BoolConstant< IsStrictlyLower_v<MT1> && IsStrictlyLower_v<MT2> && IsStrictlyLower_v<MT3> >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  YIELDSUPPER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, typename MT3 >
struct YieldsUpper<MAC,MT1,MT2,MT3>
   : public BoolConstant< IsUpper_v<MT1> && IsUpper_v<MT2> && IsUpper_v<MT3> >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  YIELDSSTRICTLYUPPER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, typename MT3 >
struct YieldsStrictlyUpper<MAC,MT1,MT2,MT3>
   : public BoolConstant< IsStrictlyUpper_v<MT1> && IsStrictlyUpper_v<MT2> && IsStrictlyUpper_v<MT3> >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  YIELDSZERO SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2, typename T3 >
struct YieldsZero<MAC,T1,T2,T3>
   : public BoolConstant< IsZero_v<T1> && IsZero_v<T2> && IsZero_v<T3> >
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
