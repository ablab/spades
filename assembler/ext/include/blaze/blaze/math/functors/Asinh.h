//=================================================================================================
/*!
//  \file blaze/math/functors/Asinh.h
//  \brief Header file for the Asinh functor
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

#ifndef _BLAZE_MATH_FUNCTORS_ASINH_H_
#define _BLAZE_MATH_FUNCTORS_ASINH_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/math/constraints/SIMDPack.h>
#include <blaze/math/shims/Asinh.h>
#include <blaze/math/simd/Asinh.h>
#include <blaze/math/typetraits/HasSIMDAsinh.h>
#include <blaze/math/typetraits/IsHermitian.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/typetraits/IsUniform.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/typetraits/IsZero.h>
#include <blaze/math/typetraits/YieldsHermitian.h>
#include <blaze/math/typetraits/YieldsLower.h>
#include <blaze/math/typetraits/YieldsStrictlyLower.h>
#include <blaze/math/typetraits/YieldsStrictlyUpper.h>
#include <blaze/math/typetraits/YieldsSymmetric.h>
#include <blaze/math/typetraits/YieldsUniform.h>
#include <blaze/math/typetraits/YieldsUpper.h>
#include <blaze/math/typetraits/YieldsZero.h>
#include <blaze/system/HostDevice.h>
#include <blaze/system/Inline.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Generic wrapper for the asinh() function.
// \ingroup functors
*/
struct Asinh
{
   //**********************************************************************************************
   /*!\brief Returns the result of the asinh() function for the given object/value.
   //
   // \param a The given object/value.
   // \return The result of the asinh() function for the given object/value.
   */
   template< typename T >
   BLAZE_ALWAYS_INLINE BLAZE_DEVICE_CALLABLE decltype(auto) operator()( T&& a ) const
   {
      return asinh( std::forward<T>( a ) );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether SIMD is enabled for the specified data type \a T.
   //
   // \return \a true in case SIMD is enabled for the data type \a T, \a false if not.
   */
   template< typename T >
   static constexpr bool simdEnabled() { return HasSIMDAsinh_v<T>; }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the operation supports padding, i.e. whether it can deal with zeros.
   //
   // \return \a true in case padding is supported, \a false if not.
   */
   static constexpr bool paddingEnabled() { return true; }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns the result of the asinh() function for the given SIMD vector.
   //
   // \param a The given SIMD vector.
   // \return The result of the asinh() function for the given SIMD vector.
   */
   template< typename T >
   BLAZE_ALWAYS_INLINE decltype(auto) load( const T& a ) const
   {
      BLAZE_CONSTRAINT_MUST_BE_SIMD_PACK( T );
      return asinh( a );
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
template< typename T >
struct YieldsUniform<Asinh,T>
   : public IsUniform<T>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  YIELDSSYMMETRIC SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT >
struct YieldsSymmetric<Asinh,MT>
   : public IsSymmetric<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  YIELDSHERMITIAN SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT >
struct YieldsHermitian<Asinh,MT>
   : public IsHermitian<MT>
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
template< typename MT >
struct YieldsLower<Asinh,MT>
   : public IsLower<MT>
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
template< typename MT >
struct YieldsStrictlyLower<Asinh,MT>
   : public IsStrictlyLower<MT>
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
template< typename MT >
struct YieldsUpper<Asinh,MT>
   : public IsUpper<MT>
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
template< typename MT >
struct YieldsStrictlyUpper<Asinh,MT>
   : public IsStrictlyUpper<MT>
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
template< typename T >
struct YieldsZero<Asinh,T>
   : public IsZero<T>
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
