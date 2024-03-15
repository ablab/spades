//=================================================================================================
/*!
//  \file blaze/math/functors/Clamp.h
//  \brief Header file for the Clamp functor
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

#ifndef _BLAZE_MATH_FUNCTORS_CLAMP_H_
#define _BLAZE_MATH_FUNCTORS_CLAMP_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/math/constraints/SIMDPack.h>
#include <blaze/math/shims/Clamp.h>
#include <blaze/math/simd/Max.h>
#include <blaze/math/simd/Min.h>
#include <blaze/math/simd/Set.h>
#include <blaze/math/typetraits/HasSIMDMax.h>
#include <blaze/math/typetraits/HasSIMDMin.h>
#include <blaze/system/HostDevice.h>
#include <blaze/system/Inline.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Generic wrapper for the clamp() function.
// \ingroup functors
*/
struct Clamp
{
   //**********************************************************************************************
   /*!\brief Returns the result of the clamp() function for the given object/value.
   //
   // \param v The given object/value to clamp.
   // \param lo The minimum to clamp \a v.
   // \param hi The maximum to clamp \a v.
   // \return The result of the clamp() function for the given object/value.
   */
   template< typename T1, typename T2, typename T3 >
   BLAZE_ALWAYS_INLINE BLAZE_DEVICE_CALLABLE decltype(auto)
      operator()( T1&& v, T2&& lo, T3&& hi ) const
   {
      return clamp( std::forward<T1>( v ), std::forward<T2>( lo ), std::forward<T3>( hi ) );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether SIMD is enabled for the specified data type \a T.
   //
   // \return \a true in case SIMD is enabled for the data type \a T, \a false if not.
   */
   template< typename T1, typename T2, typename T3 >
   static constexpr bool simdEnabled() { return HasSIMDMax_v<T1,T2> && HasSIMDMin_v<T1,T3>; }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the operation supports padding, i.e. whether it can deal with zeros.
   //
   // \return \a true in case padding is supported, \a false if not.
   */
   static constexpr bool paddingEnabled() { return true; }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns the result of the clamp() function for the given SIMD vector.
   //
   // \param v The SIMD vector to clamp.
   // \param lo The minimum to clamp \a v.
   // \param hi The maximum to clamp \a v.
   // \return The result of the clamp() function for the given SIMD vector.
   */
   template< typename T1, typename T2, typename T3 >
   BLAZE_ALWAYS_INLINE decltype(auto) load( const T1& v, const T2& lo, const T3& hi ) const
   {
      BLAZE_CONSTRAINT_MUST_BE_SIMD_PACK( T1 );
      BLAZE_CONSTRAINT_MUST_BE_SIMD_PACK( T2 );
      BLAZE_CONSTRAINT_MUST_BE_SIMD_PACK( T3 );
      return min( max( v, lo ), hi );
   }
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
