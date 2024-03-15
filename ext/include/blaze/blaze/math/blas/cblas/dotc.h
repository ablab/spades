//=================================================================================================
/*!
//  \file blaze/math/blas/cblas/dotc.h
//  \brief Header file for the CBLAS dotc wrapper functions
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

#ifndef _BLAZE_MATH_BLAS_CBLAS_DOTC_H_
#define _BLAZE_MATH_BLAS_CBLAS_DOTC_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/blas/Types.h>
#include <blaze/system/BLAS.h>
#include <blaze/util/Complex.h>
#include <blaze/util/StaticAssert.h>


namespace blaze {

//=================================================================================================
//
//  BLAS DOT PRODUCT (DOTC)
//
//=================================================================================================

//*************************************************************************************************
/*!\name BLAS dot product functions (dotc) */
//@{
#if BLAZE_BLAS_MODE

float dotc( blas_int_t n, const float* x, blas_int_t incX, const float* y, blas_int_t incY );

double dotc( blas_int_t n, const double* x, blas_int_t incX, const double* y, blas_int_t incY );

complex<float> dotc( blas_int_t n, const complex<float>* x, blas_int_t incX,
                     const complex<float>* y, blas_int_t incY );

complex<double> dotc( blas_int_t n, const complex<double>* x, blas_int_t incX,
                      const complex<double>* y, blas_int_t incY );

#endif
//@}
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a dense vector complex conjugate dot product for single precision operands
//        (\f$ s=\vec{x}*\vec{y} \f$).
// \ingroup blas
//
// \param n The size of the two dense vectors \a x and \a y \f$[0..\infty)\f$.
// \param x Pointer to the first element of vector \a x.
// \param incX The stride within vector \a x.
// \param y Pointer to the first element of vector \a y.
// \param incY The stride within vector \a y.
// \return The result of the complex conjugate dot product computation.
//
// This function performs the dense vector complex conjugate dot product for single precision
// operands based on the cblas_sdot() function (\f$ s=\vec{x}*\vec{y} \f$).
//
// \note This function can only be used if a fitting BLAS library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline float dotc( blas_int_t n, const float* x, blas_int_t incX, const float* y, blas_int_t incY )
{
   return cblas_sdot( n, x, incX, y, incY );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a dense vector complex conjugate dot product for double precision operands
//        (\f$ s=\vec{x}*\vec{y} \f$).
// \ingroup blas
//
// \param n The size of the two dense vectors \a x and \a y \f$[0..\infty)\f$.
// \param x Pointer to the first element of vector \a x.
// \param incX The stride within vector \a x.
// \param y Pointer to the first element of vector \a y.
// \param incY The stride within vector \a y.
// \return The result of the complex conjugate dot product computation.
//
// This function performs the dense vector complex conjugate dot product for double precision
// operands based on the cblas_ddot() function (\f$ s=\vec{x}*\vec{y} \f$).
//
// \note This function can only be used if a fitting BLAS library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline double dotc( blas_int_t n, const double* x, blas_int_t incX, const double* y, blas_int_t incY )
{
   return cblas_ddot( n, x, incX, y, incY );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a dense vector complex conjugate dot product for single precision
//        complex operands (\f$ s=\vec{x}*\vec{y} \f$).
// \ingroup blas
//
// \param n The size of the two dense vectors \a x and \a y \f$[0..\infty)\f$.
// \param x Pointer to the first element of vector \a x.
// \param incX The stride within vector \a x.
// \param y Pointer to the first element of vector \a y.
// \param incY The stride within vector \a y.
// \return The result of the complex conjugate dot product computation.
//
// This function performs the dense vector complex conjugate dot product for single precision
// complex operands based on the cblas_cdotc() function (\f$ s=\vec{x}*\vec{y} \f$).
//
// \note This function can only be used if a fitting BLAS library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline complex<float> dotc( blas_int_t n, const complex<float>* x, blas_int_t incX,
                            const complex<float>* y, blas_int_t incY )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<float> ) == 2UL*sizeof( float ) );

   complex<float> tmp;

#ifdef OPENBLAS_VERSION
   cblas_cdotc_sub( n, reinterpret_cast<const float*>( x ), incX,
                    reinterpret_cast<const float*>( y ), incY,
                    reinterpret_cast<openblas_complex_float*>( &tmp ) );
#else
   cblas_cdotc_sub( n, reinterpret_cast<const float*>( x ), incX,
                    reinterpret_cast<const float*>( y ), incY, &tmp );
#endif

   return tmp;
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a dense vector complex conjugate dot product for double precision
//        complex operands (\f$ s=\vec{x}*\vec{y} \f$).
// \ingroup blas
//
// \param n The size of the two dense vectors \a x and \a y \f$[0..\infty)\f$.
// \param x Pointer to the first element of vector \a x.
// \param incX The stride within vector \a x.
// \param y Pointer to the first element of vector \a y.
// \param incY The stride within vector \a y.
// \return The result of the complex conjugate dot product computation.
//
// This function performs the dense vector complex conjugate dot product for double precision
// complex operands based on the cblas_zdotc() function (\f$ s=\vec{x}*\vec{y} \f$).
//
// \note This function can only be used if a fitting BLAS library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline complex<double> dotc( blas_int_t n, const complex<double>* x, blas_int_t incX,
                             const complex<double>* y, blas_int_t incY )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<double> ) == 2UL*sizeof( double ) );

   complex<double> tmp;

#ifdef OPENBLAS_VERSION
   cblas_zdotc_sub( n, reinterpret_cast<const double*>( x ), incX,
                    reinterpret_cast<const double*>( y ), incY,
                    reinterpret_cast<openblas_complex_double*>( &tmp ) );
#else
   cblas_zdotc_sub( n, reinterpret_cast<const double*>( x ), incX,
                    reinterpret_cast<const double*>( y ), incY, &tmp );
#endif

   return tmp;
}
#endif
//*************************************************************************************************

} // namespace blaze

#endif
