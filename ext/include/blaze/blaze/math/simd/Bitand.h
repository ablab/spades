//=================================================================================================
/*!
//  \file blaze/math/simd/Bitand.h
//  \brief Header file for the SIMD bitand functionality
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

#ifndef _BLAZE_MATH_SIMD_BITAND_H_
#define _BLAZE_MATH_SIMD_BITAND_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/simd/BasicTypes.h>
#include <blaze/system/Inline.h>
#include <blaze/system/Vectorization.h>


namespace blaze {

//=================================================================================================
//
//  8-BIT INTEGRAL SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Bitwise AND ('&') of two vectors of 8-bit integral SIMD values of the same type.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return The result of the bitwise AND.
//
// This operation is only available for SSE2, AVX2, MIC and AVX512F.
*/
template< typename T >  // Type of both operands
BLAZE_ALWAYS_INLINE const T
   operator&( const SIMDi8<T>& a, const SIMDi8<T>& b ) noexcept
#if BLAZE_AVX512F_MODE || BLAZE_MIC_MODE
{
   return _mm512_and_si512( (*a).value, (*b).value );
}
#elif BLAZE_AVX2_MODE
{
   return _mm256_and_si256( (*a).value, (*b).value );
}
#elif BLAZE_SSE2_MODE
{
   return _mm_and_si128( (*a).value, (*b).value );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise AND ('&') of two vectors of 8-bit integral SIMD values of different type.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return The result of the bitwise AND.
//
// This operation is only available for SSE2, AVX2, MIC and AVX512F.
*/
template< typename T1    // Type of the left-hand side operand
        , typename T2 >  // Type of the right-hand side operand
BLAZE_ALWAYS_INLINE const SIMDuint8
   operator&( const SIMDi8<T1>& a, const SIMDi8<T2>& b ) noexcept
#if BLAZE_AVX512F_MODE || BLAZE_MIC_MODE
{
   return _mm512_and_si512( (*a).value, (*b).value );
}
#elif BLAZE_AVX2_MODE
{
   return _mm256_and_si256( (*a).value, (*b).value );
}
#elif BLAZE_SSE2_MODE
{
   return _mm_and_si128( (*a).value, (*b).value );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise AND ('&') of two vectors of 8-bit integral complex SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return The result of the bitwise AND.
//
// This operation is only available for SSE2, AVX2, MIC and AVX512F.
*/
template< typename T >  // Type of both operands
BLAZE_ALWAYS_INLINE const T
   operator&( const SIMDci8<T>& a, const SIMDci8<T>& b ) noexcept
#if BLAZE_AVX512F_MODE || BLAZE_MIC_MODE
{
   return _mm512_and_si512( (*a).value, (*b).value );
}
#elif BLAZE_AVX2_MODE
{
   return _mm256_and_si256( (*a).value, (*b).value );
}
#elif BLAZE_SSE2_MODE
{
   return _mm_and_si128( (*a).value, (*b).value );
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
/*!\brief Bitwise AND ('&') of two vectors of 16-bit integral SIMD values of the same type.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return The result of the bitwise AND.
//
// This operation is only available for SSE2, AVX2, MIC and AVX512F.
*/
template< typename T >  // Type of both operands
BLAZE_ALWAYS_INLINE const T
   operator&( const SIMDi16<T>& a, const SIMDi16<T>& b ) noexcept
#if BLAZE_AVX512F_MODE || BLAZE_MIC_MODE
{
   return _mm512_and_si512( (*a).value, (*b).value );
}
#elif BLAZE_AVX2_MODE
{
   return _mm256_and_si256( (*a).value, (*b).value );
}
#elif BLAZE_SSE2_MODE
{
   return _mm_and_si128( (*a).value, (*b).value );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise AND ('&') of two vectors of 16-bit integral SIMD values of different type.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return The result of the bitwise AND.
//
// This operation is only available for SSE2, AVX2, MIC and AVX512F.
*/
template< typename T1    // Type of the left-hand side operand
        , typename T2 >  // Type of the right-hand side operand
BLAZE_ALWAYS_INLINE const SIMDuint16
   operator&( const SIMDi16<T1>& a, const SIMDi16<T2>& b ) noexcept
#if BLAZE_AVX512F_MODE || BLAZE_MIC_MODE
{
   return _mm512_and_si512( (*a).value, (*b).value );
}
#elif BLAZE_AVX2_MODE
{
   return _mm256_and_si256( (*a).value, (*b).value );
}
#elif BLAZE_SSE2_MODE
{
   return _mm_and_si128( (*a).value, (*b).value );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise AND ('&') of two vectors of 16-bit integral complex SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return The result of the bitwise AND.
//
// This operation is only available for SSE2, AVX2, MIC and AVX512F.
*/
template< typename T >  // Type of both operands
BLAZE_ALWAYS_INLINE const T
   operator&( const SIMDci16<T>& a, const SIMDci16<T>& b ) noexcept
#if BLAZE_AVX512F_MODE || BLAZE_MIC_MODE
{
   return _mm512_and_si512( (*a).value, (*b).value );
}
#elif BLAZE_AVX2_MODE
{
   return _mm256_and_si256( (*a).value, (*b).value );
}
#elif BLAZE_SSE2_MODE
{
   return _mm_and_si128( (*a).value, (*b).value );
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
/*!\brief Bitwise AND ('&') of two vectors of 32-bit integral SIMD values of the same type.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return The result of the bitwise AND.
//
// This operation is only available for SSE2, AVX2, MIC, and AVX-512.
*/
template< typename T >  // Type of both operands
BLAZE_ALWAYS_INLINE const T
   operator&( const SIMDi32<T>& a, const SIMDi32<T>& b ) noexcept
#if BLAZE_AVX512F_MODE || BLAZE_MIC_MODE
{
   return _mm512_and_si512( (*a).value, (*b).value );
}
#elif BLAZE_AVX2_MODE
{
   return _mm256_and_si256( (*a).value, (*b).value );
}
#elif BLAZE_SSE2_MODE
{
   return _mm_and_si128( (*a).value, (*b).value );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise AND ('&') of two vectors of 32-bit integral SIMD values of different type.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return The result of the bitwise AND.
//
// This operation is only available for SSE2, AVX2, MIC, and AVX-512.
*/
template< typename T1    // Type of the left-hand side operand
        , typename T2 >  // Type of the right-hand side operand
BLAZE_ALWAYS_INLINE const SIMDuint32
   operator&( const SIMDi32<T1>& a, const SIMDi32<T2>& b ) noexcept
#if BLAZE_AVX512F_MODE || BLAZE_MIC_MODE
{
   return _mm512_and_si512( (*a).value, (*b).value );
}
#elif BLAZE_AVX2_MODE
{
   return _mm256_and_si256( (*a).value, (*b).value );
}
#elif BLAZE_SSE2_MODE
{
   return _mm_and_si128( (*a).value, (*b).value );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise AND ('&') of two vectors of 32-bit integral complex SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return The result of the bitwise AND.
//
// This operation is only available for SSE2, AVX2, MIC, and AVX-512.
*/
template< typename T >  // Type of both operands
BLAZE_ALWAYS_INLINE const T
   operator&( const SIMDci32<T>& a, const SIMDci32<T>& b ) noexcept
#if BLAZE_AVX512F_MODE || BLAZE_MIC_MODE
{
   return _mm512_and_si512( (*a).value, (*b).value );
}
#elif BLAZE_AVX2_MODE
{
   return _mm256_and_si256( (*a).value, (*b).value );
}
#elif BLAZE_SSE2_MODE
{
   return _mm_and_si128( (*a).value, (*b).value );
}
#else
= delete;
#endif
//*************************************************************************************************




//=================================================================================================
//
//  64-BIT INTEGRAL SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Bitwise AND ('&') of two vectors of 64-bit integral SIMD values of the same type.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return The result of the bitwise AND.
//
// This operation is only available for SSE2, AVX2, MIC, and AVX-512.
*/
template< typename T >  // Type of the left-hand side operand
BLAZE_ALWAYS_INLINE const T
   operator&( const SIMDi64<T>& a, const SIMDi64<T>& b ) noexcept
#if BLAZE_AVX512F_MODE || BLAZE_MIC_MODE
{
   return _mm512_and_si512( (*a).value, (*b).value );
}
#elif BLAZE_AVX2_MODE
{
   return _mm256_and_si256( (*a).value, (*b).value );
}
#elif BLAZE_SSE2_MODE
{
   return _mm_and_si128( (*a).value, (*b).value );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise AND ('&') of two vectors of 64-bit integral SIMD values of different type.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return The result of the bitwise AND.
//
// This operation is only available for SSE2, AVX2, MIC, and AVX-512.
*/
template< typename T1    // Type of the left-hand side operand
        , typename T2 >  // Type of the right-hand side operand
BLAZE_ALWAYS_INLINE const SIMDuint64
   operator&( const SIMDi64<T1>& a, const SIMDi64<T2>& b ) noexcept
#if BLAZE_AVX512F_MODE || BLAZE_MIC_MODE
{
   return _mm512_and_si512( (*a).value, (*b).value );
}
#elif BLAZE_AVX2_MODE
{
   return _mm256_and_si256( (*a).value, (*b).value );
}
#elif BLAZE_SSE2_MODE
{
   return _mm_and_si128( (*a).value, (*b).value );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise AND ('&') of two vectors of 64-bit integral complex SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return The result of the bitwise AND.
//
// This operation is only available for SSE2, AVX2, MIC, and AVX-512.
*/
template< typename T >  // Type of both operands
BLAZE_ALWAYS_INLINE const T
   operator&( const SIMDci64<T>& a, const SIMDci64<T>& b ) noexcept
#if BLAZE_AVX512F_MODE || BLAZE_MIC_MODE
{
   return _mm512_and_si512( (*a).value, (*b).value );
}
#elif BLAZE_AVX2_MODE
{
   return _mm256_and_si256( (*a).value, (*b).value );
}
#elif BLAZE_SSE2_MODE
{
   return _mm_and_si128( (*a).value, (*b).value );
}
#else
= delete;
#endif
//*************************************************************************************************

} // namespace blaze

#endif
