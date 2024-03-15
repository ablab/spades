//=================================================================================================
/*!
//  \file blaze/math/functors/Bitxor.h
//  \brief Header file for the Bitxor functor
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

#ifndef _BLAZE_MATH_FUNCTORS_BITXOR_H_
#define _BLAZE_MATH_FUNCTORS_BITXOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/math/constraints/SIMDPack.h>
#include <blaze/math/simd/Bitxor.h>
#include <blaze/math/typetraits/HasSIMDBitxor.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/typetraits/IsUniform.h>
#include <blaze/math/typetraits/YieldsSymmetric.h>
#include <blaze/math/typetraits/YieldsUniform.h>
#include <blaze/system/Inline.h>
#include <blaze/util/IntegralConstant.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Generic wrapper for the bitwise XOR ('^') operator.
// \ingroup functors
*/
struct Bitxor
{
   //**********************************************************************************************
   /*!\brief Returns the result of the bitwise XOR operator for the given objects/values.
   //
   // \param a The left-hand side object/value.
   // \param b The right-hand side object/value.
   // \return The result of the bitwise XOR operator for the given objects/values.
   */
   template< typename T1, typename T2 >
   BLAZE_ALWAYS_INLINE decltype(auto) operator()( T1&& a, T2&& b ) const
   {
      return std::forward<T1>( a ) ^ std::forward<T2>( b );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether SIMD is enabled for the specified data types \a T1 and \a T2.
   //
   // \return \a true in case SIMD is enabled for the data types \a T1 and \a T2, \a false if not.
   */
   template< typename T1, typename T2 >
   static constexpr bool simdEnabled() { return HasSIMDBitxor_v<T1,T2>; }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the operation supports padding, i.e. whether it can deal with zeros.
   //
   // \return \a true in case padding is supported, \a false if not.
   */
   static constexpr bool paddingEnabled() { return true; }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns the result of the bitwise XOR operation for the given SIMD vectors.
   //
   // \param a The left-hand side SIMD vector.
   // \param b The right-hand side SIMD vector.
   // \return The result of the bitwise XOR operation for the given SIMD vectors.
   */
   template< typename T1, typename T2 >
   BLAZE_ALWAYS_INLINE decltype(auto) load( const T1& a, const T2& b ) const
   {
      BLAZE_CONSTRAINT_MUST_BE_SIMD_PACK( T1 );
      BLAZE_CONSTRAINT_MUST_BE_SIMD_PACK( T2 );
      return a ^ b;
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
struct YieldsUniform<Bitxor,T1,T2>
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
struct YieldsSymmetric<Bitxor,MT1,MT2>
   : public BoolConstant< IsSymmetric_v<MT1> && IsSymmetric_v<MT2> >
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
