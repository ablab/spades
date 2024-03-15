//=================================================================================================
/*!
//  \file blaze/math/functors/Schur.h
//  \brief Header file for the Schur functor
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

#ifndef _BLAZE_MATH_FUNCTORS_SCHUR_H_
#define _BLAZE_MATH_FUNCTORS_SCHUR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/typetraits/IsUniform.h>
#include <blaze/math/typetraits/IsUniLower.h>
#include <blaze/math/typetraits/IsUniUpper.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/typetraits/IsZero.h>
#include <blaze/math/typetraits/YieldsLower.h>
#include <blaze/math/typetraits/YieldsStrictlyLower.h>
#include <blaze/math/typetraits/YieldsStrictlyUpper.h>
#include <blaze/math/typetraits/YieldsSymmetric.h>
#include <blaze/math/typetraits/YieldsUniform.h>
#include <blaze/math/typetraits/YieldsUniLower.h>
#include <blaze/math/typetraits/YieldsUniUpper.h>
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
/*!\brief Generic wrapper for the Schur product.
// \ingroup functors
*/
struct Schur
{
   //**********************************************************************************************
   /*!\brief Returns the result of the Schur product for the given objects/values.
   //
   // \param a The left-hand side object/value.
   // \param b The right-hand side object/value.
   // \return The result of the Schur product for the given objects/values.
   */
   template< typename T1, typename T2 >
   BLAZE_ALWAYS_INLINE BLAZE_DEVICE_CALLABLE decltype(auto) operator()( T1&& a, T2&& b ) const
   {
      return std::forward<T1>( a ) % std::forward<T2>( b );
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
template< typename T1, typename T2 >
struct YieldsUniform<Schur,T1,T2>
   : public BoolConstant< IsUniform_v<T1> && IsUniform_v<T2> >
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
template< typename MT1, typename MT2 >
struct YieldsSymmetric<Schur,MT1,MT2>
   : public BoolConstant< IsSymmetric_v<MT1> && IsSymmetric_v<MT2> >
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
template< typename MT1, typename MT2 >
struct YieldsLower<Schur,MT1,MT2>
   : public BoolConstant< IsLower_v<MT1> || IsLower_v<MT2> >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  YIELDSUNILOWER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2 >
struct YieldsUniLower<Schur,MT1,MT2>
   : public BoolConstant< IsUniLower_v<MT1> && IsUniLower_v<MT2> >
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
template< typename MT1, typename MT2 >
struct YieldsStrictlyLower<Schur,MT1,MT2>
   : public BoolConstant< IsStrictlyLower_v<MT1> || IsStrictlyLower_v<MT2> >
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
template< typename MT1, typename MT2 >
struct YieldsUpper<Schur,MT1,MT2>
   : public BoolConstant< IsUpper_v<MT1> || IsUpper_v<MT2> >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  YIELDSUNIUPPER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2 >
struct YieldsUniUpper<Schur,MT1,MT2>
   : public BoolConstant< IsUniUpper_v<MT1> && IsUniUpper_v<MT2> >
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
template< typename MT1, typename MT2 >
struct YieldsStrictlyUpper<Schur,MT1,MT2>
   : public BoolConstant< IsStrictlyUpper_v<MT1> || IsStrictlyUpper_v<MT2> >
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
template< typename T1, typename T2 >
struct YieldsZero<Schur,T1,T2>
   : public BoolConstant< IsZero_v<T1> || IsZero_v<T2> >
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
