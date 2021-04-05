//=================================================================================================
/*!
//  \file blaze/math/simd/Equal.h
//  \brief Header file for the SIMD equality functionality
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

#ifndef _BLAZE_MATH_SIMD_EQUAL_H_
#define _BLAZE_MATH_SIMD_EQUAL_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Accuracy.h>
#include <blaze/math/RelaxationFlag.h>
#include <blaze/math/simd/BasicTypes.h>
#include <blaze/system/Compiler.h>
#include <blaze/system/Inline.h>
#include <blaze/system/Vectorization.h>


namespace blaze {

//=================================================================================================
//
//  8-BIT INTEGRAL SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Equality check of two vectors of 8-bit integral SIMD values of the same type.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return \a true if the two vectors are equal, \a false if not.
//
// This operation is only available for SSE2, AVX2, MIC, and AVX-512.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename T >       // Type of both operands
BLAZE_ALWAYS_INLINE bool equal( const SIMDi8<T>& a, const SIMDi8<T>& b ) noexcept
#if BLAZE_AVX512BW_MODE
{
   return _mm512_cmpeq_epi8_mask( (*a).value, (*b).value ) == 0xffffffffffffffff;
}
#elif BLAZE_AVX2_MODE
{
   return _mm256_movemask_epi8( _mm256_cmpeq_epi8( (*a).value, (*b).value ) ) == int(0xffffffff);
}
#elif BLAZE_SSE2_MODE
{
   return _mm_movemask_epi8( _mm_cmpeq_epi8( (*a).value, (*b).value ) ) == int(0xffff);
}
#else
= delete;
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Equality check of two vectors of 8-bit integral complex SIMD values of the same type.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return \a true if the two vectors are equal, \a false if not.
//
// This operation is only available for SSE2, AVX2, MIC, and AVX-512.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename T >       // Type of both operands
BLAZE_ALWAYS_INLINE bool equal( const SIMDci8<T>& a, const SIMDci8<T>& b ) noexcept
#if BLAZE_AVX512BW_MODE
{
   return _mm512_cmpeq_epi8_mask( (*a).value, (*b).value ) == 0xffffffffffffffff;
}
#elif BLAZE_AVX2_MODE
{
   return _mm256_movemask_epi8( _mm256_cmpeq_epi8( (*a).value, (*b).value ) ) == int(0xffffffff);
}
#elif BLAZE_SSE2_MODE
{
   return _mm_movemask_epi8( _mm_cmpeq_epi8( (*a).value, (*b).value ) ) == int(0xffff);
}
#else
= delete;
#endif
/*! \endcond */
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Equality comparison of two vectors of 8-bit integral SIMD values of the same type.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return \a true if the two vectors are equal, \a false if not.
//
// This operation is only available for SSE2, AVX2, MIC, and AVX-512.
*/
template< typename T >  // Type of both operands
BLAZE_ALWAYS_INLINE bool operator==( const SIMDi8<T>& a, const SIMDi8<T>& b ) noexcept
{
   return equal<strict>( *a, *b );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison of two vectors of 8-bit integral SIMD complex values of the same type.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return \a true if the two vectors are equal, \a false if not.
//
// This operation is only available for SSE2, AVX2, MIC, and AVX-512.
*/
template< typename T >  // Type of both operands
BLAZE_ALWAYS_INLINE bool operator==( const SIMDci8<T>& a, const SIMDci8<T>& b ) noexcept
{
   return equal<strict>( *a, *b );
}
//*************************************************************************************************




//=================================================================================================
//
//  16-BIT INTEGRAL SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Equality check of two vectors of 16-bit integral SIMD values of the same type.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return \a true if the two vectors are equal, \a false if not.
//
// This operation is only available for SSE2, AVX2, MIC, and AVX-512.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename T >       // Type of both operands
BLAZE_ALWAYS_INLINE bool equal( const SIMDi16<T>& a, const SIMDi16<T>& b ) noexcept
#if BLAZE_AVX512BW_MODE
{
   return _mm512_cmpeq_epi16_mask( (*a).value, (*b).value ) == 0xffffffff;
}
#elif BLAZE_AVX2_MODE
{
   return _mm256_movemask_epi8( _mm256_cmpeq_epi16( (*a).value, (*b).value ) ) == int(0xffffffff);
}
#elif BLAZE_SSE2_MODE
{
   return _mm_movemask_epi8( _mm_cmpeq_epi16( (*a).value, (*b).value ) ) == int(0xffff);
}
#else
= delete;
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Equality check of two vectors of 16-bit integral complex SIMD values of the same type.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return \a true if the two vectors are equal, \a false if not.
//
// This operation is only available for SSE2, AVX2, MIC, and AVX-512.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename T >       // Type of both operands
BLAZE_ALWAYS_INLINE bool equal( const SIMDci16<T>& a, const SIMDci16<T>& b ) noexcept
#if BLAZE_AVX512BW_MODE
{
   return _mm512_cmpeq_epi16_mask( (*a).value, (*b).value ) == 0xffffffff;
}
#elif BLAZE_AVX2_MODE
{
   return _mm256_movemask_epi8( _mm256_cmpeq_epi16( (*a).value, (*b).value ) ) == int(0xffffffff);
}
#elif BLAZE_SSE2_MODE
{
   return _mm_movemask_epi8( _mm_cmpeq_epi16( (*a).value, (*b).value ) ) == int(0xffff);
}
#else
= delete;
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison of two vectors of 16-bit integral SIMD values of the same type.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return \a true if the two vectors are equal, \a false if not.
//
// This operation is only available for SSE2, AVX2, MIC, and AVX-512.
*/
template< typename T >  // Type of both operands
BLAZE_ALWAYS_INLINE bool operator==( const SIMDi16<T>& a, const SIMDi16<T>& b ) noexcept
{
   return equal<strict>( *a, *b );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison of two vectors of 16-bit integral complex SIMD values of the same type.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return \a true if the two vectors are equal, \a false if not.
//
// This operation is only available for SSE2, AVX2, MIC, and AVX-512.
*/
template< typename T >  // Type of both operands
BLAZE_ALWAYS_INLINE bool operator==( const SIMDci16<T>& a, const SIMDci16<T>& b ) noexcept
{
   return equal<strict>( *a, *b );
}
//*************************************************************************************************




//=================================================================================================
//
//  32-BIT INTEGRAL SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Equality check of two vectors of 32-bit integral SIMD values of the same type.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return \a true if the two vectors are equal, \a false if not.
//
// This operation is only available for SSE2, AVX2, MIC, and AVX-512.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename T >       // Type of both operands
BLAZE_ALWAYS_INLINE bool equal( const SIMDi32<T>& a, const SIMDi32<T>& b ) noexcept
#if BLAZE_AVX512F_MODE || BLAZE_MIC_MODE
{
   return _mm512_cmpeq_epi32_mask( (*a).value, (*b).value ) == 0xffff;
}
#elif BLAZE_AVX2_MODE
{
   return _mm256_movemask_epi8( _mm256_cmpeq_epi32( (*a).value, (*b).value ) ) == int(0xffffffff);
}
#elif BLAZE_SSE2_MODE
{
   return _mm_movemask_epi8( _mm_cmpeq_epi32( (*a).value, (*b).value ) ) == int(0xffff);
}
#else
= delete;
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Equality check of two vectors of 32-bit integral complex SIMD values of the same type.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return \a true if the two vectors are equal, \a false if not.
//
// This operation is only available for SSE2, AVX2, MIC, and AVX-512.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename T >       // Type of both operands
BLAZE_ALWAYS_INLINE bool equal( const SIMDci32<T>& a, const SIMDci32<T>& b ) noexcept
#if BLAZE_AVX512F_MODE || BLAZE_MIC_MODE
{
   return _mm512_cmpeq_epi32_mask( (*a).value, (*b).value ) == 0xffff;
}
#elif BLAZE_AVX2_MODE
{
   return _mm256_movemask_epi8( _mm256_cmpeq_epi32( (*a).value, (*b).value ) ) == int(0xffffffff);
}
#elif BLAZE_SSE2_MODE
{
   return _mm_movemask_epi8( _mm_cmpeq_epi32( (*a).value, (*b).value ) ) == int(0xffff);
}
#else
= delete;
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison of two vectors of 32-bit integral SIMD values of the same type.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return \a true if the two vectors are equal, \a false if not.
//
// This operation is only available for SSE2, AVX2, MIC, and AVX-512.
*/
template< typename T >  // Type of both operands
BLAZE_ALWAYS_INLINE bool operator==( const SIMDi32<T>& a, const SIMDi32<T>& b ) noexcept
{
   return equal<strict>( *a, *b );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison of two vectors of 32-bit integral complex SIMD values of the same type.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return \a true if the two vectors are equal, \a false if not.
//
// This operation is only available for SSE2, AVX2, MIC, and AVX-512.
*/
template< typename T >  // Type of both operands
BLAZE_ALWAYS_INLINE bool operator==( const SIMDci32<T>& a, const SIMDci32<T>& b ) noexcept
{
   return equal<strict>( *a, *b );
}
//*************************************************************************************************





//=================================================================================================
//
//  64-BIT INTEGRAL SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Equality check of two vectors of 64-bit integral SIMD values of the same type.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return \a true if the two vectors are equal, \a false if not.
//
// This operation is only available for SSE4.1, AVX2, MIC, and AVX-512.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename T >       // Type of both operands
BLAZE_ALWAYS_INLINE bool equal( const SIMDi64<T>& a, const SIMDi64<T>& b ) noexcept
#if BLAZE_AVX512F_MODE
{
   return _mm512_cmpeq_epi64_mask( (*a).value, (*b).value ) == 0xff;
}
#elif BLAZE_AVX2_MODE
{
   return _mm256_movemask_epi8( _mm256_cmpeq_epi64( (*a).value, (*b).value ) ) == int(0xffffffff);
}
#elif BLAZE_SSE4_MODE
{
   return _mm_movemask_epi8( _mm_cmpeq_epi64( (*a).value, (*b).value ) ) == int(0xffff);
}
#else
= delete;
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Equality check of two vectors of 64-bit integral complex SIMD values of the same type.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return \a true if the two vectors are equal, \a false if not.
//
// This operation is only available for SSE4.1, AVX2, MIC, and AVX-512.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename T >       // Type of both operands
BLAZE_ALWAYS_INLINE bool equal( const SIMDci64<T>& a, const SIMDci64<T>& b ) noexcept
#if BLAZE_AVX512F_MODE
{
   return _mm512_cmpeq_epi64_mask( (*a).value, (*b).value ) == 0xff;
}
#elif BLAZE_AVX2_MODE
{
   return _mm256_movemask_epi8( _mm256_cmpeq_epi64( (*a).value, (*b).value ) ) == int(0xffffffff);
}
#elif BLAZE_SSE4_MODE
{
   return _mm_movemask_epi8( _mm_cmpeq_epi64( (*a).value, (*b).value ) ) == int(0xffff);
}
#else
= delete;
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison of two vectors of 64-bit integral SIMD values of the same type.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return \a true if the two vectors are equal, \a false if not.
//
// This operation is only available for SSE4.1, AVX2, MIC, and AVX-512.
*/
template< typename T >  // Type of both operands
BLAZE_ALWAYS_INLINE bool operator==( const SIMDi64<T>& a, const SIMDi64<T>& b ) noexcept
{
   return equal<strict>( *a, *b );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison of two vectors of 64-bit integral SIMD complex values of the same type.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return \a true if the two vectors are equal, \a false if not.
//
// This operation is only available for SSE4.1, AVX2, MIC, and AVX-512.
*/
template< typename T >  // Type of both operands
BLAZE_ALWAYS_INLINE bool operator==( const SIMDci64<T>& a, const SIMDci64<T>& b ) noexcept
{
   return equal<strict>( *a, *b );
}
//*************************************************************************************************




//=================================================================================================
//
//  32-BIT FLOATING POINT SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Equality check of two vectors of single precision floating point SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return \a true if the two vectors are equal, \a false if not.
//
// Equal function for the comparison of two single precision floating point SIMD values. Due
// to the limited machine accuracy, a direct comparison of two floating point numbers should
// be avoided. This function offers the possibility to compare two floating-point values with
// a certain accuracy margin.
//
// This operation is only available for SSE, AVX, MIC, and AVX-512.
*/
template< RelaxationFlag RF >  // Relaxation flag
BLAZE_ALWAYS_INLINE bool equal( const SIMDfloat& a, const SIMDfloat& b ) noexcept
#if ( BLAZE_AVX512F_MODE || BLAZE_MIC_MODE ) && BLAZE_GNU_COMPILER
= delete;
#elif BLAZE_AVX512F_MODE || BLAZE_MIC_MODE
{
   if( RF == relaxed ) {
      const __m512 accu( _mm512_set1_ps( static_cast<float>( accuracy ) ) );

      const __m512 xmm1( _mm512_abs_ps( _mm512_sub_ps( a.value, b.value ) ) );
      const __m512 xmm2( _mm512_max_ps( accu, _mm512_mul_ps( accu, _mm512_abs_ps( a.value ) ) ) );
      return _mm512_cmple_ps_mask( xmm1, xmm2 ) == 0xffff;
   }
   else {
      return _mm512_cmpeq_ps_mask( a.value, b.value ) == 0xffff;
   }
}
#elif BLAZE_AVX_MODE
{
   if( RF == relaxed ) {
      const __m256 accu( _mm256_set1_ps( static_cast<float>( accuracy ) ) );
      const __m256 mask( _mm256_castsi256_ps( _mm256_set1_epi32( 0x80000000 ) ) );

      const __m256 xmm1( _mm256_andnot_ps( mask, _mm256_sub_ps( a.value, b.value ) ) );
      const __m256 xmm2( _mm256_max_ps( accu, _mm256_mul_ps( accu, _mm256_andnot_ps( mask, a.value ) ) ) );
      return _mm256_movemask_ps( _mm256_cmp_ps( xmm1, xmm2, _CMP_LE_OQ ) ) == 0xff;
   }
   else {
      return _mm256_movemask_ps( _mm256_cmp_ps( a.value, b.value, _CMP_EQ_OQ ) ) == 0xff;
   }
}
#elif BLAZE_SSE_MODE
{
   if( RF == relaxed ) {
      const __m128 accu( _mm_set1_ps( static_cast<float>( accuracy ) ) );
      const __m128 mask( _mm_castsi128_ps( _mm_set1_epi32( 0x80000000 ) ) );

      const __m128 xmm1( _mm_andnot_ps( mask, _mm_sub_ps( a.value, b.value ) ) );
      const __m128 xmm2( _mm_max_ps( accu, _mm_mul_ps( accu, _mm_andnot_ps( mask, a.value ) ) ) );
      return _mm_movemask_ps( _mm_cmple_ps( xmm1, xmm2 ) ) == 0xf;
   }
   else {
      return _mm_movemask_ps( _mm_cmpeq_ps( a.value, b.value ) ) == 0xf;
   }
}
#else
= delete;
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Equality check of two vectors of single precision complex SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return \a true if the two vectors are equal, \a false if not.
//
// Equal function for the comparison of two single precision complex SIMD values. Due to the
// limited machine accuracy, a direct comparison of two floating point numbers should be avoided.
// This function offers the possibility to compare two floating-point values with a certain
// accuracy margin.
//
// This operation is only available for SSE, AVX, MIC, and AVX-512.
*/
template< RelaxationFlag RF >  // Relaxation flag
BLAZE_ALWAYS_INLINE bool equal( const SIMDcfloat& a, const SIMDcfloat& b ) noexcept
#if ( BLAZE_AVX512F_MODE || BLAZE_MIC_MODE ) && BLAZE_GNU_COMPILER
= delete;
#elif BLAZE_AVX512F_MODE || BLAZE_MIC_MODE
{
   if( RF == relaxed ) {
      const __m512 accu( _mm512_set1_ps( static_cast<float>( accuracy ) ) );

      const __m512 xmm1( _mm512_abs_ps( _mm512_sub_ps( a.value, b.value ) ) );
      const __m512 xmm2( _mm512_max_ps( accu, _mm512_mul_ps( accu, _mm512_abs_ps( a.value ) ) ) );
      return _mm512_cmple_ps_mask( xmm1, xmm2 ) == 0xffff;
   }
   else {
      return _mm512_cmpeq_ps_mask( a.value, b.value ) == 0xffff;
   }
}
#elif BLAZE_AVX_MODE
{
   if( RF == relaxed ) {
      const __m256 accu( _mm256_set1_ps( static_cast<float>( accuracy ) ) );
      const __m256 mask( _mm256_castsi256_ps( _mm256_set1_epi32( 0x80000000 ) ) );

      const __m256 xmm1( _mm256_andnot_ps( mask, _mm256_sub_ps( a.value, b.value ) ) );
      const __m256 xmm2( _mm256_max_ps( accu, _mm256_mul_ps( accu, _mm256_andnot_ps( mask, a.value ) ) ) );
      return _mm256_movemask_ps( _mm256_cmp_ps( xmm1, xmm2, _CMP_LE_OQ ) ) == 0xff;
   }
   else {
      return _mm256_movemask_ps( _mm256_cmp_ps( a.value, b.value, _CMP_EQ_OQ ) ) == 0xff;
   }
}
#elif BLAZE_SSE_MODE
{
   if( RF == relaxed ) {
      const __m128 accu( _mm_set1_ps( static_cast<float>( accuracy ) ) );
      const __m128 mask( _mm_castsi128_ps( _mm_set1_epi32( 0x80000000 ) ) );

      const __m128 xmm1( _mm_andnot_ps( mask, _mm_sub_ps( a.value, b.value ) ) );
      const __m128 xmm2( _mm_max_ps( accu, _mm_mul_ps( accu, _mm_andnot_ps( mask, a.value ) ) ) );
      return _mm_movemask_ps( _mm_cmple_ps( xmm1, xmm2 ) ) == 0xf;
   }
   else {
      return _mm_movemask_ps( _mm_cmpeq_ps( a.value, b.value ) ) == 0xf;
   }
}
#else
= delete;
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison of two vectors of single precision floating point SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return \a true if the two vectors are equal, \a false if not.
//
// This operation is only available for SSE, AVX, MIC, and AVX-512.
*/
template< typename T1    // Type of the left-hand side operand
        , typename T2 >  // Type of the right-hand side operand
BLAZE_ALWAYS_INLINE bool operator==( const SIMDf32<T1>& a, const SIMDf32<T2>& b ) noexcept
#if ( BLAZE_AVX512F_MODE || BLAZE_MIC_MODE ) && BLAZE_GNU_COMPILER
= delete;
#elif BLAZE_AVX512F_MODE || BLAZE_MIC_MODE
{
   return _mm512_cmpeq_ps_mask( (*a).eval().value, (*b).eval().value ) == 0xffff;
}
#elif BLAZE_AVX_MODE
{
   return _mm256_movemask_ps( _mm256_cmp_ps( (*a).eval().value, (*b).eval().value, _CMP_EQ_OQ ) ) == 0xff;
}
#elif BLAZE_SSE_MODE
{
   return _mm_movemask_ps( _mm_cmpeq_ps( (*a).eval().value, (*b).eval().value ) ) == 0xf;
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison of two vectors of single precision complex SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return \a true if the two vectors are equal, \a false if not.
//
// This operation is only available for SSE, AVX, MIC, and AVX-512.
*/
BLAZE_ALWAYS_INLINE bool operator==( const SIMDcfloat& a, const SIMDcfloat& b ) noexcept
#if ( BLAZE_AVX512F_MODE || BLAZE_MIC_MODE ) && BLAZE_GNU_COMPILER
= delete;
#elif BLAZE_AVX512F_MODE || BLAZE_MIC_MODE
{
   return _mm512_cmpeq_ps_mask( a.value, b.value ) == 0xffff;
}
#elif BLAZE_AVX_MODE
{
   return _mm256_movemask_ps( _mm256_cmp_ps( a.value, b.value, _CMP_EQ_OQ ) ) == 0xff;
}
#elif BLAZE_SSE_MODE
{
   return _mm_movemask_ps( _mm_cmpeq_ps( a.value, b.value ) ) == 0xf;
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
/*! \cond BLAZE_INTERNAL */
/*!\brief Equality check of two vectors of double precision floating point SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return \a true if the two vectors are equal, \a false if not.
//
// Equal function for the comparison of two double precision floating point SIMD values. Due
// to the limited machine accuracy, a direct comparison of two floating point numbers should
// be avoided. This function offers the possibility to compare two floating-point values with
// a certain accuracy margin.
//
// This operation is only available for SSE2, AVX, MIC, and AVX-512.
*/
template< RelaxationFlag RF >  // Relaxation flag
BLAZE_ALWAYS_INLINE bool equal( const SIMDdouble& a, const SIMDdouble& b ) noexcept
#if ( BLAZE_AVX512F_MODE || BLAZE_MIC_MODE ) && BLAZE_GNU_COMPILER
= delete;
#elif BLAZE_AVX512F_MODE || BLAZE_MIC_MODE
{
   if( RF == relaxed ) {
      const __m512 accu( _mm512_set1_pd( static_cast<double>( accuracy ) ) );

      const __m512 xmm1( _mm512_abs_pd( _mm512_sub_pd( a.value, b.value ) ) );
      const __m512 xmm2( _mm512_max_pd( accu, _mm512_mul_pd( accu, _mm512_abs_pd( a.value ) ) ) );
      return _mm512_cmple_pd_mask( xmm1, xmm2 ) == 0xff;
   }
   else {
      return _mm512_cmpeq_pd_mask( a.value, b.value ) == 0xff;
   }
}
#elif BLAZE_AVX_MODE
{
   if( RF == relaxed ) {
      const __m256d accu( _mm256_set1_pd( static_cast<double>( accuracy ) ) );
      const __m256d mask( _mm256_castsi256_pd(
         _mm256_set_epi32( 0x80000000, 0x0, 0x80000000, 0x0, 0x80000000, 0x0, 0x80000000, 0x0 ) ) );

      const __m256d xmm1( _mm256_andnot_pd( mask, _mm256_sub_pd( a.value, b.value ) ) );
      const __m256d xmm2( _mm256_max_pd( accu, _mm256_mul_pd( accu, _mm256_andnot_pd( mask, a.value ) ) ) );
      return _mm256_movemask_pd( _mm256_cmp_pd( xmm1, xmm2, _CMP_LE_OQ ) ) == 0xf;
   }
   else {
      return _mm256_movemask_pd( _mm256_cmp_pd( a.value, b.value, _CMP_EQ_OQ ) ) == 0xf;
   }
}
#elif BLAZE_SSE2_MODE
{
   if( RF == relaxed ) {
      const __m128d accu( _mm_set1_pd( static_cast<double>( accuracy ) ) );
      const __m128d mask( _mm_castsi128_pd( _mm_set_epi32( 0x80000000, 0x0, 0x80000000, 0x0 ) ) );

      const __m128d xmm1( _mm_andnot_pd( mask, _mm_sub_pd( a.value, b.value ) ) );
      const __m128d xmm2( _mm_max_pd( accu, _mm_mul_pd( accu, _mm_andnot_pd( mask, a.value ) ) ) );
      return _mm_movemask_pd( _mm_cmple_pd( xmm1, xmm2 ) ) == 0x3;
   }
   else {
      return _mm_movemask_pd( _mm_cmpeq_pd( a.value, b.value ) ) == 0x3;
   }
}
#else
= delete;
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Equality check of two vectors of double precision complex SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return \a true if the two vectors are equal, \a false if not.
//
// Equal function for the comparison of two double precision complex SIMD values. Due to the
// limited machine accuracy, a direct comparison of two floating point numbers should be avoided.
// This function offers the possibility to compare two floating-point values with a certain
// accuracy margin.
//
// This operation is only available for SSE2, AVX, MIC, and AVX-512.
*/
template< RelaxationFlag RF >  // Relaxation flag
BLAZE_ALWAYS_INLINE bool equal( const SIMDcdouble& a, const SIMDcdouble& b ) noexcept
#if ( BLAZE_AVX512F_MODE || BLAZE_MIC_MODE ) && BLAZE_GNU_COMPILER
= delete;
#elif BLAZE_AVX512F_MODE || BLAZE_MIC_MODE
{
   if( RF == relaxed ) {
      const __m512 accu( _mm512_set1_pd( static_cast<double>( accuracy ) ) );

      const __m512 xmm1( _mm512_abs_pd( _mm512_sub_pd( a.value, b.value ) ) );
      const __m512 xmm2( _mm512_max_pd( accu, _mm512_mul_pd( accu, _mm512_abs_pd( a.value ) ) ) );
      return _mm512_cmple_pd_mask( xmm1, xmm2 ) == 0xff;
   }
   else {
      return _mm512_cmpeq_pd_mask( a.value, b.value ) == 0xff;
   }
}
#elif BLAZE_AVX_MODE
{
   if( RF == relaxed ) {
      const __m256d accu( _mm256_set1_pd( static_cast<double>( accuracy ) ) );
      const __m256d mask( _mm256_castsi256_pd(
         _mm256_set_epi32( 0x80000000, 0x0, 0x80000000, 0x0, 0x80000000, 0x0, 0x80000000, 0x0 ) ) );

      const __m256d xmm1( _mm256_andnot_pd( mask, _mm256_sub_pd( a.value, b.value ) ) );
      const __m256d xmm2( _mm256_max_pd( accu, _mm256_mul_pd( accu, _mm256_andnot_pd( mask, a.value ) ) ) );
      return _mm256_movemask_pd( _mm256_cmp_pd( xmm1, xmm2, _CMP_LE_OQ ) ) == 0xf;
   }
   else {
      return _mm256_movemask_pd( _mm256_cmp_pd( a.value, b.value, _CMP_EQ_OQ ) ) == 0xf;
   }
}
#elif BLAZE_SSE2_MODE
{
   if( RF == relaxed ) {
      const __m128d accu( _mm_set1_pd( static_cast<double>( accuracy ) ) );
      const __m128d mask( _mm_castsi128_pd( _mm_set_epi32( 0x80000000, 0x0, 0x80000000, 0x0 ) ) );

      const __m128d xmm1( _mm_andnot_pd( mask, _mm_sub_pd( a.value, b.value ) ) );
      const __m128d xmm2( _mm_max_pd( accu, _mm_mul_pd( accu, _mm_andnot_pd( mask, a.value ) ) ) );
      return _mm_movemask_pd( _mm_cmple_pd( xmm1, xmm2 ) ) == 0x3;
   }
   else {
      return _mm_movemask_pd( _mm_cmpeq_pd( a.value, b.value ) ) == 0x3;
   }
}
#else
= delete;
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison of two vectors of double precision floating point SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return \a true if the two vectors are equal, \a false if not.
//
// This operation is only available for SSE2, AVX, MIC, and AVX-512.
*/
template< typename T1    // Type of the left-hand side operand
        , typename T2 >  // Type of the right-hand side operand
BLAZE_ALWAYS_INLINE bool operator==( const SIMDf64<T1>& a, const SIMDf64<T2>& b ) noexcept
#if ( BLAZE_AVX512F_MODE || BLAZE_MIC_MODE ) && BLAZE_GNU_COMPILER
= delete;
#elif BLAZE_AVX512F_MODE || BLAZE_MIC_MODE
{
   return _mm512_cmpeq_pd_mask( (*a).eval().value, (*b).eval().value ) == 0xff;
}
#elif BLAZE_AVX_MODE
{
   return _mm256_movemask_pd( _mm256_cmp_pd( (*a).eval().value, (*b).eval().value, _CMP_EQ_OQ ) ) == 0xf;
}
#elif BLAZE_SSE2_MODE
{
   return _mm_movemask_pd( _mm_cmpeq_pd( (*a).eval().value, (*b).eval().value ) ) == 0x3;
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison of two vectors of double precision complex SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return \a true if the two vectors are equal, \a false if not.
//
// This operation is only available for SSE2, AVX, MIC, and AVX-512.
*/
BLAZE_ALWAYS_INLINE bool operator==( const SIMDcdouble& a, const SIMDcdouble& b ) noexcept
#if ( BLAZE_AVX512F_MODE || BLAZE_MIC_MODE ) && BLAZE_GNU_COMPILER
= delete;
#elif BLAZE_AVX512F_MODE || BLAZE_MIC_MODE
{
   return _mm512_cmpeq_pd_mask( a.value, b.value ) == 0xff;
}
#elif BLAZE_AVX_MODE
{
   return _mm256_movemask_pd( _mm256_cmp_pd( a.value, b.value, _CMP_EQ_OQ ) ) == 0xf;
}
#elif BLAZE_SSE2_MODE
{
   return _mm_movemask_pd( _mm_cmpeq_pd( a.value, b.value ) ) == 0x3;
}
#else
= delete;
#endif
//*************************************************************************************************




//=================================================================================================
//
//  INEQUALITY COMPARISON
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Inequality comparison of two SIMD vectors of the same type.
// \ingroup simd
//
// \param a The left-hand side SIMD operand.
// \param b The right-hand side SIMD operand.
// \return \a false if the two vectors are equal, \a true if not.
*/
template< typename T >  // Type of both operands
BLAZE_ALWAYS_INLINE bool operator!=( const SIMDPack<T>& a, const SIMDPack<T>& b ) noexcept
{
   return !( (*a) == (*b) );
}
//*************************************************************************************************

} // namespace blaze

#endif
