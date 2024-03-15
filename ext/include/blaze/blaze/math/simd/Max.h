//=================================================================================================
/*!
//  \file blaze/math/simd/Max.h
//  \brief Header file for the SIMD max functionality
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

#ifndef _BLAZE_MATH_SIMD_MAX_H_
#define _BLAZE_MATH_SIMD_MAX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/simd/BasicTypes.h>
#include <blaze/math/typetraits/IsSIMDPack.h>
#include <blaze/system/Inline.h>
#include <blaze/system/Vectorization.h>
#include <blaze/util/IntegralConstant.h>


namespace blaze {

//=================================================================================================
//
//  8-BIT INTEGRAL SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Componentwise maximum of two vectors of 8-bit signed integral SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return The resulting vector.
//
// This operation is only available for SSE4, AVX2, and AVX-512.
*/
BLAZE_ALWAYS_INLINE const SIMDint8 max( const SIMDint8& a, const SIMDint8& b ) noexcept
#if BLAZE_AVX512BW_MODE
{
   return _mm512_max_epi8( (*a).value, (*b).value );
}
#elif BLAZE_AVX2_MODE
{
   return _mm256_max_epi8( (*a).value, (*b).value );
}
#elif BLAZE_SSE4_MODE
{
   return _mm_max_epi8( (*a).value, (*b).value );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Componentwise maximum of two vectors of 8-bit unsigned integral SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return The resulting vector.
//
// This operation is only available for SSE2, AVX2, and AVX-512.
*/
BLAZE_ALWAYS_INLINE const SIMDuint8 max( const SIMDuint8& a, const SIMDuint8& b ) noexcept
#if BLAZE_AVX512BW_MODE
{
   return _mm512_max_epu8( (*a).value, (*b).value );
}
#elif BLAZE_AVX2_MODE
{
   return _mm256_max_epu8( (*a).value, (*b).value );
}
#elif BLAZE_SSE2_MODE
{
   return _mm_max_epu8( (*a).value, (*b).value );
}
#else
= delete;
#endif
//*************************************************************************************************




//=================================================================================================
//
//  16-BIT INTEGRAL SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Componentwise maximum of two vectors of 16-bit signed integral SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return The resulting vector.
//
// This operation is only available for SSE2, AVX2, and AVX-512.
*/
BLAZE_ALWAYS_INLINE const SIMDint16 max( const SIMDint16& a, const SIMDint16& b ) noexcept
#if BLAZE_AVX512BW_MODE
{
   return _mm512_max_epi16( (*a).value, (*b).value );
}
#elif BLAZE_AVX2_MODE
{
   return _mm256_max_epi16( (*a).value, (*b).value );
}
#elif BLAZE_SSE2_MODE
{
   return _mm_max_epi16( (*a).value, (*b).value );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Componentwise maximum of two vectors of 16-bit unsigned integral SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return The resulting vector.
//
// This operation is only available for SSE4, AVX2, and AVX-512.
*/
BLAZE_ALWAYS_INLINE const SIMDuint16 max( const SIMDuint16& a, const SIMDuint16& b ) noexcept
#if BLAZE_AVX512BW_MODE
{
   return _mm512_max_epu16( (*a).value, (*b).value );
}
#elif BLAZE_AVX2_MODE
{
   return _mm256_max_epu16( (*a).value, (*b).value );
}
#elif BLAZE_SSE4_MODE
{
   return _mm_max_epu16( (*a).value, (*b).value );
}
#else
= delete;
#endif
//*************************************************************************************************




//=================================================================================================
//
//  32-BIT INTEGRAL SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Componentwise maximim of two vectors of 32-bit signed integral SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return The resulting vector.
//
// This operation is only available for SSE4, AVX2, MIC, and AVX-512.
*/
BLAZE_ALWAYS_INLINE const SIMDint32 max( const SIMDint32& a, const SIMDint32& b ) noexcept
#if BLAZE_AVX512F_MODE || BLAZE_MIC_MODE
{
   return _mm512_max_epi32( (*a).value, (*b).value );
}
#elif BLAZE_AVX2_MODE
{
   return _mm256_max_epi32( (*a).value, (*b).value );
}
#elif BLAZE_SSE4_MODE
{
   return _mm_max_epi32( (*a).value, (*b).value );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Componentwise maximum of two vectors of 32-bit unsigned integral SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return The resulting vector.
//
// This operation is only available for SSE4, AVX2, MIC, and AVX-512.
*/
BLAZE_ALWAYS_INLINE const SIMDuint32 max( const SIMDuint32& a, const SIMDuint32& b ) noexcept
#if BLAZE_AVX512F_MODE || BLAZE_MIC_MODE
{
   return _mm512_max_epu32( (*a).value, (*b).value );
}
#elif BLAZE_AVX2_MODE
{
   return _mm256_max_epu32( (*a).value, (*b).value );
}
#elif BLAZE_SSE4_MODE
{
   return _mm_max_epu32( (*a).value, (*b).value );
}
#else
= delete;
#endif
//*************************************************************************************************




//=================================================================================================
//
//  32-BIT FLOATING POINT SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Componentwise maximum of two vectors of single precision floating point SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return The resulting vector.
//
// This operation is only available for SSE, AVX, MIC, and AVX-512.
*/
template< typename T1    // Type of the left-hand side operand
        , typename T2 >  // Type of the right-hand side operand
BLAZE_ALWAYS_INLINE const SIMDfloat
   max( const SIMDf32<T1>& a, const SIMDf32<T2>& b ) noexcept
#if BLAZE_AVX512F_MODE || BLAZE_MIC_MODE
{
   return _mm512_max_ps( (*a).eval().value, (*b).eval().value );
}
#elif BLAZE_AVX_MODE
{
   return _mm256_max_ps( (*a).eval().value, (*b).eval().value );
}
#elif BLAZE_SSE_MODE
{
   return _mm_max_ps( (*a).eval().value, (*b).eval().value );
}
#else
= delete;
#endif
//*************************************************************************************************




//=================================================================================================
//
//  64-BIT FLOATING POINT SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Componentwise maximum of two vectors of double precision floating point SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return The resulting vector.
//
// This operation is only available for SSE2, AVX, MIC, and AVX-512.
*/
template< typename T1    // Type of the left-hand side operand
        , typename T2 >  // Type of the right-hand side operand
BLAZE_ALWAYS_INLINE const SIMDdouble
   max( const SIMDf64<T1>& a, const SIMDf64<T2>& b ) noexcept
#if BLAZE_AVX512F_MODE || BLAZE_MIC_MODE
{
   return _mm512_max_pd( (*a).eval().value, (*b).eval().value );
}
#elif BLAZE_AVX_MODE
{
   return _mm256_max_pd( (*a).eval().value, (*b).eval().value );
}
#elif BLAZE_SSE2_MODE
{
   return _mm_max_pd( (*a).eval().value, (*b).eval().value );
}
#else
= delete;
#endif
//*************************************************************************************************

} // namespace blaze

#endif
